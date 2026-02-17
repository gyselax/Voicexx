// SPDX-License-Identifier: MIT
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "chargedensitycalculator.hpp"
#include "collisions_inter.hpp"
#include "collisions_intra.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "fem_1d_poisson_solver.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry_xvx.hpp"
#include "input.hpp"
#include "krook_source_adaptive.hpp"
#include "krook_source_constant.hpp"
#include "maxwellian_source.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "qnsolver.hpp"
#include "restartinitialisation.hpp"
#include "sheath.yaml.hpp"
#include "singlemodeperturbinitialisation.hpp"
#include "species_info.hpp"
#include "species_init.hpp"
#include "spline_definitions_xvx.hpp"
#include "spline_interpolator.hpp"
#include "splitrighthandsidesolver.hpp"
#include "splitvlasovsolver.hpp"

using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

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
        std::cout << "SETTING FILE!" << std::endl;
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
            SplineInterpPointsX>(conf_voicexx, "x");
    IdxRangeVx const mesh_vx = init_spline_dependent_idx_range<
            GridVx,
            BSplinesVx,
            SplineInterpPointsVx>(conf_voicexx, "vx");
    IdxRangeXVx const meshXVx(mesh_x, mesh_vx);

    IdxRangeSp const idx_range_kinsp = init_species(conf_voicexx);

    IdxRangeSpXVx const meshSpXVx(idx_range_kinsp, meshXVx);
    IdxRangeSpVx const meshSpVx(idx_range_kinsp, mesh_vx);

    SplineXBuilder const builder_x(mesh_x);
    SplineVxBuilder const builder_vx(mesh_vx);

    // Initialisation of the distribution function
    DFieldMemSpVx allfequilibrium(meshSpVx);
    MaxwellianEquilibrium const init_fequilibrium
            = maxwellian_equilibrium::init_from_input(idx_range_kinsp, conf_voicexx);
    init_fequilibrium(get_field(allfequilibrium));

    ddc::expose_to_pdi("iter_start", iter_start);

    DFieldMemSpXVx allfdistribu(meshSpXVx);
    double time_start(0);
    if (iter_start == 0) {
        SingleModePerturbInitialisation const init = SingleModePerturbInitialisation::
                init_from_input(get_const_field(allfequilibrium), idx_range_kinsp, conf_voicexx);
        init(get_field(allfdistribu));
    } else {
        RestartInitialisation const restart(time_start);
        restart(get_field(allfdistribu));
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
            throw std::invalid_argument("Invalid krook name, allowed values are: 'constant', or "
                                        "'adaptive'.");
        }
    }

    // Maxwellian source
    MaxwellianSource const rhs_plasma_source
            = maxwellian_source::init_from_input(meshSpXVx, conf_voicexx);
    rhs_operators.emplace_back(rhs_plasma_source);

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

    DFieldMemVx const quadrature_coeffs(neumann_spline_quadrature_coefficients<
                                        Kokkos::DefaultExecutionSpace>(mesh_vx, builder_vx));

    ChargeDensityCalculator rhs(get_const_field(quadrature_coeffs));

    // Create the objects needed for the Poisson solver. These objects must not go out of scope
    // until after the simulation has run.
#if defined(PERIODIC_RDIMX) && !defined(INPUT_MESH)
    FFTPoissonSolver<IdxRangeX, IdxRangeX, Kokkos::DefaultExecutionSpace> poisson_solver(mesh_x);
#else
    FEM1DPoissonSolver poisson_solver(builder_x, spline_x_evaluator);
#endif
    QNSolver const poisson(poisson_solver, rhs);

    PredCorr const predcorr(boltzmann, poisson);

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
    ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium_host);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(get_field(allfdistribu), time_start, deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
