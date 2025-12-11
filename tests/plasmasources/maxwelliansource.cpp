// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <pdi.h>

#include "geometry.hpp"
#include "maxwellian_source.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"
#include "trapezoid_quadrature.hpp"

TEST(MaxwellianSource, Moments)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_size(100);

    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IdxStepVx const vx_size(30);

    IdxStepSp const nb_species(2);
    IdxRangeSp const idx_range_sp(IdxSp(0), nb_species);
    IdxSp const my_iion = idx_range_sp.front();
    IdxSp const my_ielec = idx_range_sp.back();

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX idx_range_x(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx idx_range_vx(SplineInterpPointsVx::get_domain<GridVx>());

    SplineXBuilder const builder_x(idx_range_x);
    SplineVxBuilder const builder_vx(idx_range_vx);

    host_t<DFieldMemX> quadrature_coeffs_x
            = trapezoid_quadrature_coefficients<Kokkos::DefaultHostExecutionSpace>(idx_range_x);
    host_t<DFieldMemVx> quadrature_coeffs_vx
            = trapezoid_quadrature_coefficients<Kokkos::DefaultHostExecutionSpace>(idx_range_vx);
    host_t<Quadrature<IdxRangeX>> const integrate_x(get_const_field(quadrature_coeffs_x));
    host_t<Quadrature<IdxRangeVx>> const integrate_v(get_const_field(quadrature_coeffs_vx));

    {
        host_t<DFieldMemSp> charges(idx_range_sp);
        charges(my_ielec) = -1.;
        charges(my_iion) = 1.;
        host_t<DFieldMemSp> masses(idx_range_sp);
        ddc::parallel_fill(masses, 1.);

        // Initialisation of the distribution function
        ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));
    }

    DFieldMemSpXVx allfdistribu(IdxRangeSpXVx(idx_range_sp, idx_range_x, idx_range_vx));

    // Initialisation of the distribution function
    ddc::parallel_fill(allfdistribu, 0.);

    // Maxwellian source test
    double const extent_source = 0.2;
    double const stiffness_source = 0.1;
    double const amplitude_source = 1.;
    double const density_source = 1;
    double const temperature_source_elec = 0.5;
    double const temperature_source_ions = 2.0;
    double const deltat = 1.;

    MaxwellianSource const maxw_source(
            idx_range_x,
            idx_range_vx,
            idx_range_sp,
            extent_source,
            stiffness_source,
            amplitude_source,
            density_source,
            temperature_source_elec,
            temperature_source_ions);

    maxw_source(get_field(allfdistribu), deltat);

    IdxRangeSpX idx_range_spx(idx_range_sp, idx_range_x);
    host_t<DFieldMemSpX> density(idx_range_spx);
    host_t<DFieldMemSpX> fluid_velocity(idx_range_spx);
    host_t<DFieldMemSpX> temperature(idx_range_spx);

    host_t<DFieldMemVx> values_density(idx_range_vx);
    host_t<DFieldMemVx> values_fluid_velocity(idx_range_vx);
    host_t<DFieldMemVx> values_temperature(idx_range_vx);
    ddc::for_each(idx_range_spx, [&](IdxSpX const ispx) {
        // density
        ddc::parallel_deepcopy(values_density, allfdistribu[ispx]);
        density(ispx)
                = integrate_v(Kokkos::DefaultHostExecutionSpace(), get_const_field(values_density));

        // fluid velocity
        ddc::for_each(idx_range_vx, [&](IdxVx const iv) {
            values_fluid_velocity(iv) = values_density(iv) * ddc::coordinate(iv);
        });
        fluid_velocity(ispx) = integrate_v(
                                       Kokkos::DefaultHostExecutionSpace(),
                                       get_const_field(values_fluid_velocity))
                               / density(ispx);

        // temperature
        ddc::for_each(idx_range_vx, [&](IdxVx const iv) {
            values_temperature(iv)
                    = values_density(iv) * std::pow(ddc::coordinate(iv) - fluid_velocity(ispx), 2);
        });
        temperature(ispx) = integrate_v(
                                    Kokkos::DefaultHostExecutionSpace(),
                                    get_const_field(values_temperature))
                            / density(ispx);
    });

    // source amplitude
    double error_source_amplitude_ions
            = integrate_x(Kokkos::DefaultHostExecutionSpace(), get_const_field(density[my_iion]))
              - amplitude_source;
    double error_source_amplitude_elec
            = integrate_x(Kokkos::DefaultHostExecutionSpace(), get_const_field(density[my_ielec]))
              - amplitude_source;
    EXPECT_LE(error_source_amplitude_ions, 1e-8);
    EXPECT_LE(error_source_amplitude_elec, 1e-8);

    double error_fluid_velocity_elec(0);
    double error_temperature_elec(0);
    ddc::for_each(idx_range_x, [&](IdxX const ix) {
        error_fluid_velocity_elec
                = std::fmax(std::fabs(fluid_velocity(my_ielec, ix)), error_fluid_velocity_elec);
        error_temperature_elec = std::
                fmax(std::fabs(temperature(my_ielec, ix) - temperature_source_elec),
                     error_temperature_elec);
    });
    EXPECT_LE(error_fluid_velocity_elec, 1e-8);
    EXPECT_LE(error_temperature_elec, 1e-8);

    double error_fluid_velocity_ions(0);
    double error_temperature_ions(0);
    ddc::for_each(idx_range_x, [&](IdxX const ix) {
        error_fluid_velocity_ions
                = std::fmax(std::fabs(fluid_velocity(my_iion, ix)), error_fluid_velocity_ions);
        error_temperature_ions = std::
                fmax(std::fabs(temperature(my_iion, ix) - temperature_source_ions),
                     error_temperature_ions);
    });
    EXPECT_LE(error_fluid_velocity_ions, 1e-8);
    EXPECT_LE(
            error_temperature_ions,
            1e-3); // the tolerance is lower because of the higher temperature

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
