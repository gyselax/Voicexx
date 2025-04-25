// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "charge_exchange.hpp"
#include "constantfluidinitialisation.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "geometry_moments.hpp"
#include "ionisation.hpp"
#include "recombination.hpp"
#include "species_info.hpp"


/**
 * This test initialises a plasma and a fluid species with flat profiles,
 * compute the reaction rates and compare them to an analytical solution.
*/
static void TestDiffusiveNeutralsRateCoefficients()
{
    CoordX const x_min(0.0);
    CoordX const x_max(10);
    IdxStepX const x_size(512);

    CoordVx const vx_min(-8);
    CoordVx const vx_max(8);
    IdxStepVx const vx_size(50);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX meshX(SplineInterpPointsX::get_domain<GridX>());

    // Kinetic and neutral species index range initialisation
    IdxStepSp const nb_kinspecies(2);
    IdxRangeSp const idx_range_kinsp(IdxSp(0), nb_kinspecies);

    IdxStepSp const nb_fluidspecies(1);
    IdxRangeSp const idx_range_fluidsp(IdxSp(idx_range_kinsp.back() + 1), nb_fluidspecies);

    IdxRangeSp const idx_range_allsp(IdxSp(0), nb_kinspecies + nb_fluidspecies);

    host_t<DFieldMemSp> masses(idx_range_allsp);
    host_t<DFieldSp> kinetic_masses = masses[idx_range_kinsp];
    host_t<DFieldSp> fluid_masses = masses[idx_range_fluidsp];

    host_t<DFieldMemSp> charges(idx_range_allsp);
    host_t<DFieldSp> kinetic_charges = charges[idx_range_kinsp];
    host_t<DFieldSp> fluid_charges = charges[idx_range_fluidsp];

    IdxSp const my_iion = idx_range_kinsp.front();
    IdxSp const my_ielec = idx_range_kinsp.back();
    IdxSp const my_ifluid = idx_range_fluidsp.front();

    kinetic_charges(my_ielec) = -1.;
    kinetic_charges(my_iion) = 1.;

    double const mass_ion(400.), mass_elec(1.);
    kinetic_masses(my_ielec) = mass_elec;
    kinetic_masses(my_iion) = mass_ion;

    // neutrals charge is zero
    fluid_charges(my_ifluid) = 0.;

    double const neutral_mass(1.);
    fluid_masses(my_ifluid) = neutral_mass;

    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    // Moments index range initialisation
    IdxStepMom const nb_fluid_moments(1);
    IdxRangeMom const meshM(GeometryMX::density_idx, nb_fluid_moments);
    ddc::init_discrete_space<GridMom>();

    IdxRangeSpX idx_range_fluidspx = IdxRangeSpX(idx_range_fluidsp, meshX);

    {
        // we verify that we recover the right k_cx_0
        double const n_0 = 1e20;
        double const T_0_test = 2.2704067636837;
        double const precision = 1e-13;
        ChargeExchangeRate const charge_exchange(n_0, T_0_test);
        EXPECT_NEAR(charge_exchange.get_Kcx0(), 1., precision);
    }


    ChargeExchangeRate charge_exchange(1.);
    IonisationRate ionisation(1.);
    RecombinationRate recombination(1.);

    DFieldMemSpMomX neutrals_alloc(IdxRangeSpMomX(idx_range_fluidsp, meshM, meshX));
    DFieldSpMomX neutrals = get_field(neutrals_alloc);

    host_t<DFieldMemSpMom> moments_init(IdxRangeSpMom(idx_range_fluidsp, meshM));
    ddc::parallel_fill(moments_init, 1.);
    ConstantFluidInitialisation<GridX> fluid_init(get_const_field(moments_init));
    fluid_init(neutrals);

    DFieldMemSpX kinsp_density_alloc(IdxRangeSpX(idx_range_kinsp, meshX));
    DFieldMemSpX kinsp_velocity_alloc(IdxRangeSpX(idx_range_kinsp, meshX));
    DFieldMemSpX kinsp_temperature_alloc(IdxRangeSpX(idx_range_kinsp, meshX));

    DFieldSpX kinsp_density = get_field(kinsp_density_alloc);
    DFieldSpX kinsp_velocity = get_field(kinsp_velocity_alloc);
    DFieldSpX kinsp_temperature = get_field(kinsp_temperature_alloc);

    double const kinsp_density_eq(1.);
    double const kinsp_velocity_eq(0.0);
    double const kinsp_temperature_eq(1.);
    ddc::parallel_fill(kinsp_density, kinsp_density_eq);
    ddc::parallel_fill(kinsp_velocity, kinsp_velocity_eq);
    ddc::parallel_fill(kinsp_temperature, kinsp_temperature_eq);

    // building reaction rates
    DFieldMemSpX charge_exchange_rate_alloc(idx_range_fluidspx);
    DFieldMemSpX ionisation_rate_alloc(idx_range_fluidspx);
    DFieldMemSpX recombination_rate_alloc(idx_range_fluidspx);

    DFieldSpX charge_exchange_rate = get_field(charge_exchange_rate_alloc);
    DFieldSpX ionisation_rate = get_field(ionisation_rate_alloc);
    DFieldSpX recombination_rate = get_field(recombination_rate_alloc);

    charge_exchange(
            charge_exchange_rate,
            get_const_field(kinsp_density),
            get_const_field(kinsp_temperature));
    ionisation(ionisation_rate, get_const_field(kinsp_density), get_const_field(kinsp_temperature));
    recombination(
            recombination_rate,
            get_const_field(kinsp_density),
            get_const_field(kinsp_temperature));

    double mean_cx_rate = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            idx_range_fluidspx,
            0.,
            ddc::reducer::sum<double>(),
            charge_exchange_rate);

    double mean_i_rate = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            idx_range_fluidspx,
            0.,
            ddc::reducer::sum<double>(),
            ionisation_rate);

    double mean_r_rate = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            idx_range_fluidspx,
            0.,
            ddc::reducer::sum<double>(),
            recombination_rate);

    mean_cx_rate /= meshX.size();
    mean_i_rate /= meshX.size();
    mean_r_rate /= meshX.size();

    EXPECT_NEAR(mean_cx_rate, 1.783406341061044, 1e-13);
    EXPECT_NEAR(mean_i_rate, 1.130359390036803, 1e-13);
    EXPECT_NEAR(mean_r_rate, 7.638123065868132e-06, 1e-13);

    { // we verify that the choice of normalisation does not change the result
        double const n_0 = 1e20;
        double const T_0 = 10;
        // 13 and 7 are chosen randomly. It should never change the result.
        ChargeExchangeRate const cx_rate1(n_0, T_0);
        for (double n_factor = 1e-5; n_factor < 1e6; n_factor *= 10) {
            for (double T_factor = 0.1; T_factor < 1e3; T_factor *= 10) {
                ChargeExchangeRate const cx_rate2(n_factor * n_0, T_factor * T_0);
                IdxSpX const origin(0, 0);
                IdxRangeSpX const test_range(origin, IdxStep<Species, GridX>(1, 1));
                DFieldMemSpX cx_result1_alloc(test_range);
                DFieldSpX cx_result1(get_field(cx_result1_alloc));
                DFieldMemSpX cx_result2_alloc(test_range);
                DFieldSpX cx_result2(get_field(cx_result2_alloc));
                auto result_1_host = ddc::create_mirror_view_and_copy(cx_result1);
                auto result_2_host = ddc::create_mirror_view_and_copy(cx_result2);
                cx_rate1(
                        cx_result1,
                        get_const_field(kinsp_density),
                        get_const_field(kinsp_temperature));
                cx_rate2(
                        cx_result2,
                        get_const_field(kinsp_density),
                        get_const_field(kinsp_temperature));
                EXPECT_NEAR(result_1_host(origin), result_2_host(origin), 1e-13);
            }
        }
    }

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}

TEST(ReactionRate, NeutralsRateCoefficients)
{
    TestDiffusiveNeutralsRateCoefficients();
}
