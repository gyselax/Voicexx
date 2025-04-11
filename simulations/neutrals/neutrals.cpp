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

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "charge_exchange.hpp"
#include "chargedensitycalculator.hpp"
#include "collisions_inter.hpp"
#include "collisions_intra.hpp"
#include "constantfluidinitialisation.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "fem_1d_poisson_solver.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "ionisation.hpp"
#include "irighthandside.hpp"
#include "kinetic_source.hpp"
#include "krook_source_adaptive.hpp"
#include "krook_source_constant.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "neutrals.yml.hpp"
#include "noenergytransfercoupling.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "pdi_out_neutrals.yml.hpp"
#include "predcorr_hybrid.hpp"
#include "qnsolver.hpp"
#include "recombination.hpp"
#include "restartinitialisationwithneutrals.hpp"
#include "samegridfluidsolver.hpp"
#include "singlemodeperturbinitialisation.hpp"
#include "species_info.hpp"
#include "species_init.hpp"
#include "spline_interpolator.hpp"
#include "splitrighthandsidesolver.hpp"
#include "splitvlasovsolver.hpp"

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
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
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
    // --> Mesh info
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
    SplineXBuilder const builder_x(meshXVx);
    SplineVxBuilder const builder_vx(meshXVx);
    SplineVxBuilder_1d const builder_vx_poisson(mesh_vx);

    IdxRangeSp idx_range_kinsp;
    IdxRangeSp idx_range_fluidsp;
    init_species_withfluid(idx_range_kinsp, idx_range_fluidsp, conf_voicexx);

    // Initialisation of kinetic species distribution function
    IdxRangeSpVx const meshSpVx(idx_range_kinsp, mesh_vx);
    DFieldMemSpVx allfequilibrium(meshSpVx);
    MaxwellianEquilibrium const init_fequilibrium
            = MaxwellianEquilibrium::init_from_input(idx_range_kinsp, conf_voicexx);
    init_fequilibrium(get_field(allfequilibrium));

    ddc::expose_to_pdi("iter_start", iter_start);

    double time_start(0);
    IdxRangeSpXVx const meshSpXVx(idx_range_kinsp, meshXVx);
    DFieldMemSpXVx allfdistribu(meshSpXVx);

    // Moments index range initialisation
    IdxStepMom const nb_fluid_moments(1);
    IdxRangeMom const meshM(IdxMom(0), nb_fluid_moments);
    ddc::init_discrete_space<GridMom>();
    // Neutral species initialisation
    DFieldMemSpMomX neutrals_alloc(IdxRangeSpMomX(idx_range_fluidsp, meshM, mesh_x));
    DFieldSpMomX neutrals = get_field(neutrals_alloc);

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
        ConstantFluidInitialisation fluid_init(get_const_field(moments_init_host));
        fluid_init(neutrals);

    } else {
        RestartInitialisationWithNeutrals<GridX> const restart(iter_start, time_start);
        restart(get_field(allfdistribu), get_field(neutrals));
    }
    auto allfequilibrium_host = ddc::create_mirror_view_and_copy(get_field(allfequilibrium));

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
    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);
    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

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

    DFieldMemVx const quadrature_coeffs_alloc(
            neumann_spline_quadrature_coefficients<
                    Kokkos::DefaultExecutionSpace>(mesh_vx, builder_vx_poisson));
    ChargeDensityCalculator rhs(get_const_field(quadrature_coeffs_alloc));

    // Create the objects needed for the Poisson solver. These objects must not go out of scope
    // until after the simulation has run.
#if defined(PERIODIC_RDIMX) && !defined(INPUT_MESH)
    FFTPoissonSolver<IdxRangeX, IdxRangeX, Kokkos::DefaultExecutionSpace> poisson_solver(mesh_x);
#else
    SplineXBuilder_1d const builder_x_poisson(mesh_x);
    SplineXEvaluator_1d const spline_x_evaluator_poisson(bv_x_min, bv_x_max);
    FEM1DPoissonSolver poisson_solver(builder_x_poisson, spline_x_evaluator_poisson);
#endif
    QNSolver const poisson(poisson_solver, rhs);

    // Initialisation of the neutrals
    double const normalisation_coeff
            = PCpp_double(conf_voicexx, ".DiffusiveNeutralSolver.normalisation_coeff_neutrals");
    double const norm_coeff_rate
            = PCpp_double(conf_voicexx, ".DiffusiveNeutralSolver.norm_coeff_rate_neutrals");

    // The CX coefficient needs to be first constructed in order to write a correct initstate file. Check pdi_out_neutrals.yml.hpp for a closer look.
    ChargeExchangeRate charge_exchange(norm_coeff_rate);
    IonisationRate ionisation(norm_coeff_rate);
    RecombinationRate recombination(norm_coeff_rate);

    SplineXBuilder_1d const spline_x_builder_neutrals(mesh_x);
    SplineXEvaluator_1d const spline_x_evaluator_neutrals(bv_x_min, bv_x_max);

    DFieldMemVx const quadrature_coeffs_neutrals(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(mesh_vx));

    double const neutrals_wall_extent = PCpp_double(conf_voicexx, ".NeutralKrook.extent");
    double const neutrals_wall_stiffness = PCpp_double(conf_voicexx, ".NeutralKrook.stiffness");
    double const neutrals_wall_amplitude = PCpp_double(conf_voicexx, ".NeutralKrook.amplitude");

    SameGridFluidSolver const neutralsolver(
            charge_exchange,
            ionisation,
            recombination,
            normalisation_coeff,
            spline_x_builder_neutrals,
            spline_x_evaluator_neutrals,
            get_const_field(quadrature_coeffs_neutrals),
            neutrals_wall_extent,
            neutrals_wall_stiffness,
            neutrals_wall_amplitude,
            mesh_x);

    NoEnergyExchangeCoupling const kineticfluidcoupling(
            PCpp_double(conf_voicexx, ".KineticFluidCouplingSource.density_coupling_coeff"),
            PCpp_double(conf_voicexx, ".KineticFluidCouplingSource.momentum_coupling_coeff"),
            PCpp_double(conf_voicexx, ".KineticFluidCouplingSource.energy_coupling_coeff"),
            ionisation,
            recombination,
            normalisation_coeff,
            get_const_field(quadrature_coeffs_alloc),
            neutrals_wall_extent,
            neutrals_wall_stiffness,
            mesh_x);

    PredCorrHybrid const predcorr(boltzmann, neutralsolver, poisson, kineticfluidcoupling);

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", ddc::discrete_space<BSplinesX>().ncells());
    ddc::expose_to_pdi("Nvx_spline_cells", ddc::discrete_space<BSplinesVx>().ncells());
    expose_mesh_to_pdi("MeshX", mesh_x);
    expose_mesh_to_pdi("MeshVx", mesh_vx);
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
    ddc::expose_to_pdi("normalisation_coeff_neutrals", normalisation_coeff);
    ddc::expose_to_pdi("norm_coeff_rate_neutrals", norm_coeff_rate);
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
