// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Lagrange_interpolator.hpp"
#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "chargedensitycalculator.hpp"
#include "constantfluidinitialisation.hpp"
#include "constantrate.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "fem_1d_poisson_solver.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "noenergytransfercoupling.hpp"
#include "nullfluidsolver.hpp"
#include "predcorr.hpp"
#include "predcorr_hybrid.hpp"
#include "qnsolver.hpp"
#include "singlemodeperturbinitialisation.hpp"
#include "species_info.hpp"
#include "spline_interpolator.hpp"
#include "splitvlasovsolver.hpp"

/**
 * This test creates one instance of the PredCorr class, and one instance of the 
 * PredCorrFluid class. With the NullFluidSolver class, 
 * this test verifies that the distribution function is the same after applying the 
 * predictor corrector scheme to it, with and without the fluid species.
 */
TEST(GeometryXM, PredCorrHybrid)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_ncells(50);

    CoordVx const vx_min(-8);
    CoordVx const vx_max(8);
    IdxStepVx const vx_ncells(50);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_ncells);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_ncells);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX meshX(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx meshVx(SplineInterpPointsVx::get_domain<GridVx>());
    IdxRangeXVx meshXVx(meshX, meshVx);

    SplineXBuilder const builder_x(meshX);
    SplineVxBuilder const builder_vx(meshVx);

    // Kinetic species index range initialisation
    IdxStepSp const nb_kinspecies(2);
    IdxRangeSp const idx_range_kinsp(IdxSp(0), nb_kinspecies);

    IdxSp const iion = idx_range_kinsp.front();
    IdxSp const ielec = idx_range_kinsp.back();

    host_t<DFieldMemSp> kinetic_charges(idx_range_kinsp);
    kinetic_charges(ielec) = -1.;
    kinetic_charges(iion) = 1.;

    host_t<DFieldMemSp> kinetic_masses(idx_range_kinsp);
    double const mass_ion(400.), mass_elec(1.);
    kinetic_masses(ielec) = mass_elec;
    kinetic_masses(iion) = mass_ion;

    // Fluid species index range initialisation
    IdxStepSp const nb_fluidspecies(1);
    IdxRangeSp const idx_range_fluidsp(IdxSp(idx_range_kinsp.back() + 1), nb_fluidspecies);

    // Fluid charges
    host_t<DFieldMemSp> fluid_charges(idx_range_fluidsp);
    ddc::parallel_fill(fluid_charges, 0.);

    host_t<DFieldMemSp> fluid_masses(idx_range_fluidsp);
    ddc::parallel_fill(fluid_masses, mass_ion);

    // Create the index range of kinetic species + fluid species
    IdxRangeSp const idx_range_allsp(IdxSp(0), nb_kinspecies + nb_fluidspecies);

    // Create a Field that contains charges of all species
    host_t<DFieldMemSp> charges(idx_range_allsp);

    // fill the Field with charges of kinetic species
    for (IdxSp isp : idx_range_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }

    // fill the Field with charges of fluid species
    for (IdxSp isp : idx_range_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }

    // Create a Field that contains masses of kinetic and fluid species
    host_t<DFieldMemSp> masses(idx_range_allsp);

    // fill the Field with masses of kinetic species
    for (IdxSp isp : idx_range_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }

    // fill the Field with masses of fluid species
    for (IdxSp isp : idx_range_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }

    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    // Initialisation of kinetic species distribution function
    DFieldMemSpXVx allfdistribu_alloc(IdxRangeSpXVx(idx_range_kinsp, meshX, meshVx));
    DFieldSpXVx allfdistribu = get_field(allfdistribu_alloc);

    host_t<DFieldMemSp> kinsp_density_eq(idx_range_kinsp);
    host_t<DFieldMemSp> kinsp_velocity_eq(idx_range_kinsp);
    host_t<DFieldMemSp> kinsp_temperature_eq(idx_range_kinsp);

    ddc::parallel_fill(kinsp_density_eq, 1.);
    ddc::parallel_fill(kinsp_velocity_eq, 0.);
    ddc::parallel_fill(kinsp_temperature_eq, 1.);

    DFieldMemSpVx allfequilibrium_alloc(IdxRangeSpVx(idx_range_kinsp, meshVx));
    DFieldSpVx allfequilibrium = get_field(allfequilibrium_alloc);
    MaxwellianEquilibrium const init_fequilibrium(
            std::move(kinsp_density_eq),
            std::move(kinsp_temperature_eq),
            std::move(kinsp_velocity_eq));
    init_fequilibrium(allfequilibrium);

    host_t<FieldMemSp<int>> init_perturb_mode(idx_range_kinsp);
    ddc::parallel_fill(init_perturb_mode, 2);
    host_t<DFieldMemSp> init_perturb_amplitude(idx_range_kinsp);
    ddc::parallel_fill(init_perturb_amplitude, 0.1);

    SingleModePerturbInitialisation const
            init(get_const_field(allfequilibrium),
                 std::move(init_perturb_mode),
                 std::move(init_perturb_amplitude));
    init(allfdistribu);

    // Moments index range initialisation
    IdxStepMom const nb_fluid_moments(3);
    IdxRangeMom const meshM(IdxMom(0), nb_fluid_moments);
    ddc::init_discrete_space<GridMom>();

    IdxMom idensity(0);
    IdxMom iflux(1);
    IdxMom istress(2);

    // Initialisation of fluid species moments
    DFieldMemSpMomX fluid_moments_alloc(IdxRangeSpMomX(idx_range_fluidsp, meshM, meshX));
    DFieldSpMomX fluid_moments = get_field(fluid_moments_alloc);

    host_t<DFieldMemSpMom> moments_init(IdxRangeSpMom(idx_range_fluidsp, meshM));
    ddc::parallel_fill(moments_init[idensity], 1.);
    ddc::parallel_fill(moments_init[iflux], 0.);
    ddc::parallel_fill(moments_init[istress], 1.);

    ConstantFluidInitialisation<GridX> fluid_init(get_const_field(moments_init));
    fluid_init(fluid_moments);

#ifdef PERIODIC_RDIMX
    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
#else
    ddc::ConstantExtrapolationRule<X> bv_x_min(x_min);
    ddc::ConstantExtrapolationRule<X> bv_x_max(x_max);
#endif

    // Creating operators
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
    PreallocatableSplineInterpolator const
            spline_x_interpolator(builder_x, spline_x_evaluator, meshXVx);

    IdxStepVx static constexpr gwvx {0};
    LagrangeInterpolator<GridVx, BCond::DIRICHLET, BCond::DIRICHLET, GridX, GridVx> const
            lagrange_vx_non_preallocatable_interpolator(3, gwvx);
    PreallocatableLagrangeInterpolator<
            GridVx,
            BCond::DIRICHLET,
            BCond::DIRICHLET,
            GridX,
            GridVx> const lagrange_vx_interpolator(lagrange_vx_non_preallocatable_interpolator);

    BslAdvectionSpatial<GeometryXVx, GridX> const advection_x(spline_x_interpolator);
    BslAdvectionVelocity<GeometryXVx, GridVx> const advection_vx(lagrange_vx_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    DFieldMemVx const quadrature_coeffs = neumann_spline_quadrature_coefficients<
            Kokkos::DefaultExecutionSpace>(meshVx, builder_vx);

    ChargeDensityCalculator rhs(get_const_field(quadrature_coeffs));
#ifdef PERIODIC_RDIMX
    FFTPoissonSolver<IdxRangeX> poisson_solver(meshX);
#else
    FEM1DPoissonSolver const poisson_solver(builder_x, spline_x_evaluator);
#endif
    QNSolver const poisson(poisson_solver, rhs);

    ConstantRate const charge_exchange(0.0);
    ConstantRate const ionisation(0.0);
    ConstantRate const recombination(0.0);
    double const normalisation_coeff(1.0);

    // kinetic fluid coupling term
    double const neutrals_wall_extent = 0.;
    double const neutrals_wall_stiffness = 1.;
    NoEnergyExchangeCoupling const kineticfluidcoupling(
            1.0,
            0.0,
            0.0,
            ionisation,
            recombination,
            normalisation_coeff,
            get_const_field(quadrature_coeffs),
            neutrals_wall_extent,
            neutrals_wall_stiffness,
            meshX);

    // construction of predcorr without fluid species
    PredCorr const predcorr(vlasov, poisson);

    // distribution function to be evolved by predcorr without fluid species
    DFieldMemSpXVx allfdistribu_predcorr_alloc(get_idx_range(allfdistribu));
    DFieldSpXVx allfdistribu_predcorr = get_field(allfdistribu_predcorr_alloc);
    ddc::parallel_deepcopy(allfdistribu_predcorr, allfdistribu);

    double const time_start(0.);
    int const nb_iter(10);
    double const deltat(0.1);
    predcorr(allfdistribu_predcorr, time_start, deltat, nb_iter);

    // construction of predcorr with fluid species
    NullFluidSolver<GridX> const fluidsolver(idx_range_fluidsp);
    PredCorrHybrid<GridX> const predcorr_hybrid(vlasov, fluidsolver, poisson, kineticfluidcoupling);
    predcorr_hybrid(allfdistribu, fluid_moments, time_start, deltat, nb_iter);

    auto allfdistribu_host = ddc::create_mirror_view_and_copy(allfdistribu);

    auto allfdistribu_predcorr_host = ddc::create_mirror_view_and_copy(allfdistribu_predcorr);

    /**
     * Since the fluid model uses NullFluidSolver, 
     * the distribution function evolved with PredCorrFluid and PredCorr 
     * should be equal
     */
    double const tolerance(1.e-12);
    ddc::for_each(get_idx_range(allfdistribu), [&](IdxSpXVx const ispxvx) {
        EXPECT_LE(
                std::fabs(allfdistribu_host(ispxvx) - allfdistribu_predcorr_host(ispxvx)),
                tolerance);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
