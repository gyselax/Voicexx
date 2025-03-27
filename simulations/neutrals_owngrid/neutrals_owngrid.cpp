// SPDX-License-Identifier: MIT
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>
#include <ddc/pdi.hpp>

#include <paraconf.h>

#include "ddc/discrete_space.hpp"
#include "ddc/uniform_point_sampling.hpp"

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "chargedensitycalculator.hpp"
#include "collisions_inter.hpp"
#include "collisions_intra.hpp"
#include "constantfluidinitialisation.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "densitycoupling.hpp"
#include "fem_1d_poisson_solver.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "geometry_neutrals.hpp"
#include "igridneutralcoupling.hpp"
#include "input.hpp"
#include "iplasmaneutralscoupling.hpp"
#include "irighthandside.hpp"
#include "kinetic_source.hpp"
#include "krook_source_adaptive.hpp"
#include "krook_source_constant.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "neutrals_owngrid.yml.hpp"
#include "nullfluidsolver.hpp"
#include "nullplasmaneutralcoupling.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "qnsolver.hpp"
#include "singlemodeperturbinitialisation.hpp"
#include "species_info.hpp"
#include "species_init.hpp"
#include "spline_interpolator.hpp"
#include "splitrighthandsidesolver.hpp"
#include "splitvlasovsolver.hpp"

// used for the fact that we have a different grid for the neutrals
#include "charge_exchange.hpp"
#include "densitycoupling.hpp"
#include "diffgridsfluidsolver.hpp"
#include "ionisation.hpp"
#include "nullplasmaneutralscoupling.hpp"
#include "pdi_out_neutrals_owngrid.yaml.hpp"
#include "predcorr_hybrid.hpp"
#include "recombination.hpp"
#include "restartinitialisationwithneutrals.hpp"

using std::chrono::steady_clock;

int main(int argc, char** argv)
{
    long int iter_start;
    PC_tree_t conf_voicexx;
    parse_executable_arguments(
            conf_voicexx,
            iter_start,
            argc,
            argv,
            (std::string(mesh_params_yaml) + std::string(params_yaml)).c_str());
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG); // loading pdi_out
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

    if constexpr (
            (ddcHelper::is_non_uniform_interpolation_points_v<SplineInterpPointsX>)
            || (ddc::is_non_uniform_bsplines_v<BSplinesX>)
            || (ddc::is_non_uniform_bsplines_v<BSplinesVx>)) {
        std::string spline_mesh_filename(PCpp_string(conf_voicexx, ".SplineMesh.grid_file"));
        size_t spline_mesh_filename_size = spline_mesh_filename.size();
        PDI_multi_expose(
                "setFilename",
                "filename_size",
                &spline_mesh_filename_size,
                PDI_OUT,
                "filename",
                spline_mesh_filename.c_str(),
                PDI_OUT,
                NULL);
    }

    // Reading config
    // --> mesh for the plasma
    IdxRangeX const mesh_x = init_spline_dependent_idx_range<
            GridX,
            BSplinesX,
            SplineInterpPointsX>(conf_voicexx, "x"); // constructing the mesh in X
    IdxRangeVx const mesh_vx = init_spline_dependent_idx_range<
            GridVx,
            BSplinesVx,
            SplineInterpPointsVx>(conf_voicexx, "vx"); // and the mesh in v
    IdxRangeXVx const meshXVx(mesh_x, mesh_vx); //merging the two

    // Initialisation of the spline builders
    SplineXBuilder const builder_x(mesh_x);
    SplineVxBuilder const builder_vx(mesh_vx);

    // Initialise the two IdxRangeSp for kinetic species and fluid species
    IdxRangeSp idx_range_kinsp;
    IdxRangeSp idx_range_fluidsp;
    init_species_withfluid(idx_range_kinsp, idx_range_fluidsp, conf_voicexx);

    // Initialisation of kinetic species distribution function
    IdxRangeSpVx const meshSpVx(idx_range_kinsp, mesh_vx);
    DFieldMemSpVx allfequilibrium(meshSpVx);
    MaxwellianEquilibrium const init_fequilibrium
            = MaxwellianEquilibrium::init_from_input(idx_range_kinsp, conf_voicexx);
    init_fequilibrium(get_field(allfequilibrium));
    auto allfequilibrium_host = ddc::create_mirror_view_and_copy(get_field(allfequilibrium));

    ddc::expose_to_pdi("iter_start", iter_start);

    double time_start(0);
    IdxRangeSpXVx const meshSpXVx(idx_range_kinsp, meshXVx);
    DFieldMemSpXVx allfdistribu(meshSpXVx);

    // Moments index range initialisation
    IdxStepMom const nb_fluid_moments(1);
    IdxRangeMom const meshM(IdxMom(0), nb_fluid_moments);
    ddc::init_discrete_space<GridMom>();

    // Neutral species initialisation
    // We begin by constructing the neutral mesh
    Coord<X> min(PCpp_double(conf_voicexx, ".NeutralMesh.x_min"));
    Coord<X> max(PCpp_double(conf_voicexx, ".NeutralMesh.x_max"));
    IdxStep<GridXNeutrals> ncells(PCpp_int(conf_voicexx, ".NeutralMesh.x_ncells"));
    ddc::init_discrete_space<BSplinesXNeutrals>(min, max, ncells);
    ddc::init_discrete_space<GridXNeutrals>(
            SplineInterpPointsXNeutrals::get_sampling<GridXNeutrals>());
    IdxRangeXn mesh_x_neutrals = SplineInterpPointsXNeutrals::get_domain<GridXNeutrals>();
    DFieldMemSpMomXn neutrals_alloc(IdxRangeSpMomXn(idx_range_fluidsp, meshM, mesh_x_neutrals));
    DFieldSpMomXn neutrals = get_field(neutrals_alloc);

    if (iter_start == 0) { // if we start a new simulation
        // we need to add a perturbation otherwise it will stay at equilibrium
        SingleModePerturbInitialisation const init = SingleModePerturbInitialisation::
                init_from_input(get_const_field(allfequilibrium), idx_range_kinsp, conf_voicexx);
        init(get_field(allfdistribu));

        //neutrals get init according to the input file
        host_t<DFieldMemSpMom> moments_init_host(IdxRangeSpMom(idx_range_fluidsp, meshM));
        for (IdxSp const isp : idx_range_fluidsp) {
            PC_tree_t const conf_nisp = PCpp_get(
                    conf_voicexx,
                    ".NeutralSpeciesInfo[%d]",
                    (isp - idx_range_fluidsp.front()).value());
            ddc::parallel_fill(moments_init_host[isp], PCpp_double(conf_nisp, ".density_eq"));
        }
        ConstantFluidInitialisation<GridXNeutrals> fluid_init(get_const_field(moments_init_host));
        fluid_init(neutrals);

    } else {
        RestartInitialisationWithNeutrals<GridXNeutrals> const restart(iter_start, time_start);
        restart(get_field(allfdistribu), get_field(neutrals));
    }

    // --> Algorithm info
    double const deltat = PCpp_double(conf_voicexx, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_voicexx, ".Algorithm.nbiter"));

    // --> Output info
    double const time_diag = PCpp_double(conf_voicexx, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

    // Creating operators
#ifdef PERIODIC_RDIMX
    // Macro is required to ensure spline_x_evaluator stays in the scope
    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
#else
    ddc::ConstantExtrapolationRule<X> bv_x_min(ddc::coordinate(mesh_x.front()));
    ddc::ConstantExtrapolationRule<X> bv_x_max(ddc::coordinate(mesh_x.back()));
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
#endif

    ddc::ConstantExtrapolationRule<Vx> bv_vx_min(ddc::coordinate(mesh_vx.front()));
    ddc::ConstantExtrapolationRule<Vx> bv_vx_max(ddc::coordinate(mesh_vx.back()));
    SplineVxEvaluator const spline_vx_evaluator(bv_vx_min, bv_vx_max);
    PreallocatableSplineInterpolator const
            spline_x_interpolator(builder_x, spline_x_evaluator, meshXVx);
    PreallocatableSplineInterpolator const
            spline_vx_interpolator(builder_vx, spline_vx_evaluator, meshXVx);

    BslAdvectionSpatial<GeometryXVx, GridX> const advection_x(spline_x_interpolator);
    BslAdvectionVelocity<GeometryXVx, GridVx> const advection_vx(spline_vx_interpolator);

    // list of rhs operators
    std::vector<std::reference_wrapper<IRightHandSide const>> rhs_operators;
    std::vector<KrookSourceConstant> krook_source_constant_vector;
    std::vector<KrookSourceAdaptive> krook_source_adaptive_vector;
    // Krook operators initialisation
    int const nb_rhsKrook(PCpp_len(conf_voicexx, ".Krook"));
    for (int ik = 0; ik < nb_rhsKrook; ++ik) {
        // --> Krook info
        PC_tree_t const conf_krook = PCpp_get(conf_voicexx, ".Krook[%d]", ik);

        static std::map<std::string, RhsType>
                str2rhstype {{"source", RhsType::Source}, {"sink", RhsType::Sink}};
        RhsType type = str2rhstype[PCpp_string(conf_krook, ".type")];
        std::string const krook_name = PCpp_string(conf_krook, ".name");
        if (krook_name == "constant") {
            krook_source_constant_vector.emplace_back(
                    mesh_x,
                    mesh_vx,
                    type,
                    PCpp_double(conf_krook, ".extent"),
                    PCpp_double(conf_krook, ".stiffness"),
                    PCpp_double(conf_krook, ".amplitude"),
                    PCpp_double(conf_krook, ".density"),
                    PCpp_double(conf_krook, ".temperature"));
            rhs_operators.emplace_back(krook_source_constant_vector.back());

        } else if (krook_name == "adaptive") {
            krook_source_adaptive_vector.emplace_back(
                    mesh_x,
                    mesh_vx,
                    type,
                    PCpp_double(conf_krook, ".extent"),
                    PCpp_double(conf_krook, ".stiffness"),
                    PCpp_double(conf_krook, ".amplitude"),
                    PCpp_double(conf_krook, ".density"),
                    PCpp_double(conf_krook, ".temperature"));
            rhs_operators.emplace_back(krook_source_adaptive_vector.back());
        } else {
            throw std::invalid_argument(
                    "Invalid krook name, allowed values are: 'constant', or 'adaptive'.");
        }
    }

    // Kinetic source
    KineticSource const rhs_kinetic_source(
            mesh_x,
            mesh_vx,
            PCpp_double(conf_voicexx, ".KineticSource.extent"),
            PCpp_double(conf_voicexx, ".KineticSource.stiffness"),
            PCpp_double(conf_voicexx, ".KineticSource.amplitude"),
            PCpp_double(conf_voicexx, ".KineticSource.density"),
            PCpp_double(conf_voicexx, ".KineticSource.energy"),
            PCpp_double(conf_voicexx, ".KineticSource.temperature"));
    rhs_operators.emplace_back(rhs_kinetic_source);

    CollisionsIntra const
            collisions_intra(meshSpXVx, PCpp_double(conf_voicexx, ".CollisionsInfo.nustar0"));
    rhs_operators.emplace_back(collisions_intra);

    std::optional<CollisionsInter> collisions_inter;
    if (PCpp_bool(conf_voicexx, ".CollisionsInfo.enable_inter")) {
        collisions_inter.emplace(meshSpXVx, PCpp_double(conf_voicexx, ".CollisionsInfo.nustar0"));
        rhs_operators.emplace_back(*collisions_inter);
    }
    SplitVlasovSolver const vlasov(advection_x, advection_vx);
    SplitRightHandSideSolver const boltzmann(vlasov, rhs_operators);

    DFieldMemVx const quadrature_coeffs_alloc(neumann_spline_quadrature_coefficients<
                                              Kokkos::DefaultExecutionSpace>(mesh_vx, builder_vx));
    ChargeDensityCalculator rhs(get_const_field(quadrature_coeffs_alloc));

    // Create the objects needed for the Poisson solver. These objects must not go out of scope
    // until after the simulation has run.
#if defined(PERIODIC_RDIMX) && !defined(INPUT_MESH)
    FFTPoissonSolver<IdxRangeX, IdxRangeX, Kokkos::DefaultExecutionSpace> poisson_solver(mesh_x);
#else
    FEM1DPoissonSolver poisson_solver(builder_x, spline_x_evaluator);
#endif
    QNSolver const poisson(poisson_solver, rhs);

    // Initialisation of the neutrals
    double const mean_free_path = PCpp_double(conf_voicexx, ".DiffusiveSolver.mean_free_path");
    double const temperature_normalisation = PCpp_double(conf_voicexx, ".DiffusiveSolver.T_0");
    double const density_normalisation = PCpp_double(conf_voicexx, ".DiffusiveSolver.n_0");

    // The CX coefficient needs to be first constructed in order to write a correct initstate file. Check pdi_out_neutrals.yml.hpp for a closer look.
    ChargeExchangeRate charge_exchange(density_normalisation, temperature_normalisation);
    IonisationRate ionisation(
            density_normalisation,
            temperature_normalisation,
            charge_exchange.get_Kcx0());
    RecombinationRate recombination(
            density_normalisation,
            temperature_normalisation,
            charge_exchange.get_Kcx0());
    printf(" <====================================================>\n"
           " The K_{cx,0} taken for this simulation is %fe-14.\n"
           " It has been computed from n_0=%e and T_0=%e.\n"
           " However the L_{cx,0} taken is %e.\n"
           " <====================================================>\n",
           charge_exchange.get_Kcx0(),
           density_normalisation,
           temperature_normalisation,
           mean_free_path);

    // splines to interpolate from one grid to another
    SplineXNeutralsBuilder const spline_builder_on_Xn(mesh_x_neutrals);
    SplineXn_GridXnEvaluator const spline_evaluator_on_Xn(bv_x_min, bv_x_max);
    SplineX_GridXnEvaluator const interpolator_from_X_to_Xn(bv_x_min, bv_x_max);

    DFieldMemVx const quadrature_coeffs_neutrals(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(mesh_vx));

    // depending if we want to solve the transport for the neutral species
    // we choose the corresponding neutral solver
    std::unique_ptr<IFluidSolver<GridXNeutrals>> ptr_neutral_solver;
    if (PCpp_bool(conf_voicexx, ".DiffusiveSolver.on")) {
        std::string bc_flux_input
                = PCpp_string(conf_voicexx, ".DiffusiveSolver.boundary_condition");
        DiffGridsFluidSolver::NeutralFluxBoundaryCondition bc_flux
                = DiffGridsFluidSolver::neutral_flux_boundary_condition(bc_flux_input);

        ptr_neutral_solver = std::make_unique<DiffGridsFluidSolver>(
                charge_exchange,
                ionisation,
                recombination,
                mean_free_path,
                spline_builder_on_Xn,
                spline_evaluator_on_Xn,
                builder_x,
                interpolator_from_X_to_Xn,
                get_const_field(quadrature_coeffs_neutrals),
                bc_flux,
                PCpp_double(conf_voicexx, ".DiffusiveSolver.recycling_coefficient"));
    } else {
        ptr_neutral_solver = std::make_unique<NullFluidSolver<GridXNeutrals>>(idx_range_fluidsp);
    }

    // depending if we want the plasma and the neutrals to exchange,
    // we choose the coupling
    std::unique_ptr<DensityCoupling> ptr_kinfluidcoupling;
    if (PCpp_bool(conf_voicexx, ".KineticFluidCoupling.on")) {
        ptr_kinfluidcoupling = std::make_unique<DensityCoupling>(
                PCpp_double(conf_voicexx, ".KineticFluidCoupling.density_coupling_coeff"),
                PCpp_double(conf_voicexx, ".KineticFluidCoupling.momentum_coupling_coeff"),
                PCpp_double(conf_voicexx, ".KineticFluidCoupling.energy_coupling_coeff"),
                ionisation,
                recombination,
                builder_x,
                interpolator_from_X_to_Xn,
                mean_free_path,
                get_const_field(quadrature_coeffs_alloc));
    } else {
        ptr_kinfluidcoupling = std::make_unique<NullPlasmaNeutralsCoupling<GridXNeutrals>>();
    }

    PredCorrHybrid<GridXNeutrals> const
            predcorr(boltzmann, *ptr_neutral_solver, poisson, *ptr_kinfluidcoupling);

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", ddc::discrete_space<BSplinesX>().ncells());
    ddc::expose_to_pdi("Nvx_spline_cells", ddc::discrete_space<BSplinesVx>().ncells());
    expose_mesh_to_pdi("MeshX", mesh_x);
    expose_mesh_to_pdi("MeshVx", mesh_vx);
    expose_mesh_to_pdi("MeshNeutrals", mesh_x_neutrals);
    ddc::expose_to_pdi("Lx", ddcHelper::total_interval_length(mesh_x));
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("Nkinspecies", idx_range_kinsp.size());
    ddc::expose_to_pdi(
            "fdistribu_charges",
            ddc::discrete_space<Species>().charges()[idx_range_kinsp]);
    ddc::expose_to_pdi(
            "fdistribu_masses",
            ddc::discrete_space<Species>().masses()[idx_range_kinsp]);
    ddc::expose_to_pdi(
            "neutrals_masses",
            ddc::discrete_space<Species>().masses()[idx_range_fluidsp]);
    ddc::expose_to_pdi("temperature_normalisation", temperature_normalisation);
    ddc::expose_to_pdi("density_normalisation", density_normalisation);
    ddc::expose_to_pdi("k_cx_0", charge_exchange.get_Kcx0());
    ddc::expose_to_pdi("mean_free_path", mean_free_path);
    ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium_host);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(get_field(allfdistribu), neutrals, time_start, deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
