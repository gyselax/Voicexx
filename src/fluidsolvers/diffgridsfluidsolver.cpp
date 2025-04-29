// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "diffgridsfluidsolver.hpp"
#include "geometry.hpp"
#include "geometry_neutrals.hpp"
#include "ireactionrate.hpp"
#include "rk2.hpp"
#include "species_info.hpp"

DiffGridsFluidSolver::DiffGridsFluidSolver(
        IReactionRate const& charge_exchange,
        IReactionRate const& ionisation,
        IReactionRate const& recombination,
        double const mean_free_path,
        SplineXNeutralsBuilder const& spline_builder_on_Xn,
        SplineXn_GridXnEvaluator const& spline_evaluator_on_Xn,
        SplineXBuilder const& spline_builder_on_X,
        SplineX_GridXnEvaluator const& intepolator_from_X_to_Xn,
        DConstFieldVx const& quadrature_coeffs)

    : m_charge_exchange(charge_exchange)
    , m_ionisation(ionisation)
    , m_recombination(recombination)
    , m_mean_free_path(mean_free_path)
    , m_spline_builder_on_Xn(spline_builder_on_Xn)
    , m_spline_evaluator_on_Xn(spline_evaluator_on_Xn)
    , m_spline_builder_on_X(spline_builder_on_X)
    , m_interpolator_from_X_to_Xn(intepolator_from_X_to_Xn)
    , m_quadrature_coeffs(quadrature_coeffs)
{
}

IdxSp DiffGridsFluidSolver::find_ion(IdxRangeSp const idx_range_kinsp) const
{
    bool ion_found = false;
    IdxSp iion;
    for (IdxSp const isp : idx_range_kinsp) {
        if (charge(isp) > 0.) {
            ion_found = true;
            iion = isp;
        }
    }
    if (!ion_found) {
        throw std::runtime_error("ion not found");
    }
    assert(idx_range_kinsp.size() == 2);

    return iion;
}

void DiffGridsFluidSolver::get_derivative(
        DFieldSpMomXn dn,
        DConstFieldSpMomXn neutrals,
        DConstFieldSpX density,
        DConstFieldSpX velocity,
        DConstFieldSpX temperature) const
{
    IdxRangeSpXn idx_range_fluidspxn(get_idx_range(neutrals));
    IdxRangeSpX idx_range_ratespx(get_idx_range<Species>(neutrals), get_idx_range<GridX>(density));

    // building reaction rates
    // pay attention, these are normalised to Kcx0=10^-14
    DFieldMemSpX charge_exchange_rate_alloc(idx_range_ratespx);
    DFieldMemSpX ionisation_rate_alloc(idx_range_ratespx);
    DFieldMemSpX recombination_rate_alloc(idx_range_ratespx);

    DFieldSpX charge_exchange_rate = get_field(charge_exchange_rate_alloc);
    DFieldSpX ionisation_rate = get_field(ionisation_rate_alloc);
    DFieldSpX recombination_rate = get_field(recombination_rate_alloc);

    m_charge_exchange(charge_exchange_rate, density, temperature);
    m_ionisation(ionisation_rate, density, temperature);
    m_recombination(recombination_rate, density, temperature);

    // expose to pdi the reaction coefficients
    auto cx_host = ddc::create_mirror_view_and_copy(charge_exchange_rate);
    auto i_host = ddc::create_mirror_view_and_copy(ionisation_rate);
    auto r_host = ddc::create_mirror_view_and_copy(recombination_rate);
    ddc::PdiEvent("reaction_rate_expose")
            .with("charge_exchange_rate", cx_host)
            .with("ionisation_rate", i_host)
            .with("recombination_rate", r_host);

    // create the fields to interpolate the plasma quantities on the neutral grid
    IdxRangeSpXn idx_range_kinspxn(
            get_idx_range<Species>(density),
            get_idx_range<GridXNeutrals>(neutrals));
    DFieldMemSpXn density_on_Xn_alloc(idx_range_kinspxn);
    DFieldSpXn density_on_Xn = get_field(density_on_Xn_alloc);
    DFieldMemSpXn velocity_on_Xn_alloc(idx_range_kinspxn);
    DFieldSpXn velocity_on_Xn = get_field(velocity_on_Xn_alloc);
    DFieldMemSpXn temperature_on_Xn_alloc(idx_range_kinspxn);
    DFieldSpXn temperature_on_Xn = get_field(temperature_on_Xn_alloc);

    DFieldMemSpXn cx_rate_Xn_alloc(idx_range_fluidspxn);
    DFieldMemSpXn i_rate_Xn_alloc(idx_range_fluidspxn);
    DFieldMemSpXn r_rate_Xn_alloc(idx_range_fluidspxn);
    DFieldSpXn charge_exchange_rate_on_Xn = get_field(cx_rate_Xn_alloc);
    DFieldSpXn ionisation_rate_on_Xn = get_field(i_rate_Xn_alloc);
    DFieldSpXn recombination_rate_on_Xn = get_field(r_rate_Xn_alloc);

    // do the interpolation on the neutral grid
    interpolate_on_neutral_grid(density_on_Xn, density);
    interpolate_on_neutral_grid(velocity_on_Xn, velocity);
    interpolate_on_neutral_grid(temperature_on_Xn, temperature);

    interpolate_on_neutral_grid(charge_exchange_rate_on_Xn, get_const_field(charge_exchange_rate));
    interpolate_on_neutral_grid(ionisation_rate_on_Xn, get_const_field(ionisation_rate));
    interpolate_on_neutral_grid(recombination_rate_on_Xn, get_const_field(recombination_rate));

    // compute diffusive model equation terms
    DFieldMemSpXn density_equilibrium_alloc(idx_range_fluidspxn);
    DFieldMemSpXn density_equilibrium_velocity_alloc(idx_range_fluidspxn);
    DFieldMemSpXn neutral_pressure_alloc(idx_range_fluidspxn);
    DFieldMemSpXn pressure_diffusion_coefficient_alloc(idx_range_fluidspxn);
    DFieldSpXn density_equilibrium = get_field(density_equilibrium_alloc); // n_eq
    DFieldSpXn density_equilibrium_velocity
            = get_field(density_equilibrium_velocity_alloc); // n_eq * u_i
    DFieldSpXn neutral_pressure
            = get_field(neutral_pressure_alloc); // p_N=T_i*n_N because (T_i=T_N)
    DFieldSpXn pressure_diffusion_coefficient
            = get_field(pressure_diffusion_coefficient_alloc); // \hat D

    IdxSp const iion(find_ion(get_idx_range<Species>(density)));
    IdxMom const ineutral_density(0);
    double const sqrt_mass_ratio(Kokkos::sqrt(mass(ielec()) / mass(iion)));
    double mean_free_path = m_mean_free_path;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_fluidspxn,
            KOKKOS_LAMBDA(IdxSpXn const ifspxn) {
                IdxSp const isp(ifspxn);
                IdxXn const ixn(ifspxn);

                //in this for loop we construct the terms n_eq*u_i, p_N and \hat D
                //we will take their derivatives afterwards
                double density_neutrals = neutrals(ifspxn, ineutral_density);
                double density_elec = density_on_Xn(ielec(), ixn);
                double density_ions = density_on_Xn(iion, ixn);
                double K_cx = charge_exchange_rate_on_Xn(ifspxn);
                double K_i = ionisation_rate_on_Xn(ifspxn);
                double K_r = recombination_rate_on_Xn(ifspxn);
                neutral_pressure(ifspxn) = temperature_on_Xn(iion, ixn) * density_neutrals;
                density_equilibrium(ifspxn) = (density_neutrals * density_ions * K_cx
                                               + density_ions * density_elec * K_r)
                                              / (density_ions * K_cx + density_elec * K_i);
                density_equilibrium_velocity(ifspxn)
                        = density_equilibrium(ifspxn) * velocity_on_Xn(iion, ixn);
                pressure_diffusion_coefficient(ifspxn)
                        = mean_free_path / (density_ions * K_cx + density_elec * K_i);
                // density source is not solved here, we only solve transport.
            });

    // compute the gradient of the neutral pressure
    DFieldMemSpXn gradx_pressure_alloc(idx_range_fluidspxn);
    DFieldSpXn gradx_pressure = get_field(gradx_pressure_alloc);
    ddc::for_each(get_idx_range<Species>(neutrals), [&](IdxSp const isp) {
        DBSFieldMemXn pressure_spline_x_coeff(get_spline_idx_range(m_spline_builder_on_Xn));
        m_spline_builder_on_Xn(
                get_field(pressure_spline_x_coeff),
                get_const_field(neutral_pressure[isp]));
        m_spline_evaluator_on_Xn
                .deriv(gradx_pressure[isp], get_const_field(pressure_spline_x_coeff));
    });

    // compute the neutral particle flux
    DFieldMemSpXn particle_flux_alloc(idx_range_fluidspxn);
    DFieldSpXn particle_flux = get_field(particle_flux_alloc); // \Gamma_N
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_fluidspxn,
            KOKKOS_LAMBDA(IdxSpXn const ifspxn) {
                particle_flux(ifspxn) = density_equilibrium_velocity(ifspxn) // convection
                                        - pressure_diffusion_coefficient(ifspxn)
                                                  * gradx_pressure(ifspxn); //diffusion
            });

    // enforce the boundary condition
    IdxRangeXn range_fluid_x(idx_range_fluidspxn);
    IdxXn min = range_fluid_x.front();
    IdxXn max = range_fluid_x.back();
    ddc::parallel_fill(particle_flux[min], 0);
    ddc::parallel_fill(particle_flux[max], 0);

    // compute the gradient of the flux
    // currently, splines are used but this will be changed
    DFieldMemSpXn grad_particle_flux_alloc(idx_range_fluidspxn);
    DFieldSpXn grad_particle_flux = get_field(grad_particle_flux_alloc); // \partial_x \Gamma_N
    ddc::for_each(get_idx_range<Species>(neutrals), [&](IdxSp const isp) {
        DBSFieldMemXn flux_spline_x_coeff(get_spline_idx_range(m_spline_builder_on_Xn));
        m_spline_builder_on_Xn(get_field(flux_spline_x_coeff), get_const_field(particle_flux[isp]));
        m_spline_evaluator_on_Xn
                .deriv(grad_particle_flux[isp], get_const_field(flux_spline_x_coeff));
    });

    // compute the neutral derivative
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_fluidspxn,
            KOKKOS_LAMBDA(IdxSpXn const ifspxn) {
                dn(ifspxn, ineutral_density) = -sqrt_mass_ratio * grad_particle_flux(ifspxn);
            }); // density source is not solved here, we only solve transport.

    // we expose to pdi the coefficients
    host_t<DFieldMemSpXn> n_eq_ui_host(get_idx_range(density_equilibrium_velocity));
    ddc::parallel_deepcopy(n_eq_ui_host, density_equilibrium_velocity);
    host_t<DFieldMemSpXn> pressure_grad_host(get_idx_range(gradx_pressure));
    ddc::parallel_deepcopy(pressure_grad_host, gradx_pressure);
    host_t<DFieldMemSpXn> diffusion_coeff_host(get_idx_range(pressure_diffusion_coefficient));
    ddc::parallel_deepcopy(diffusion_coeff_host, pressure_diffusion_coefficient);
    host_t<DFieldMemSpXn> particle_flux_host(get_idx_range(particle_flux));
    ddc::parallel_deepcopy(particle_flux_host, particle_flux);
    host_t<DFieldMemSpXn> flux_grad_host(get_idx_range(grad_particle_flux));
    ddc::parallel_deepcopy(flux_grad_host, grad_particle_flux);
    ddc::PdiEvent("diff_conv_expose")
            .with("pressure_grad", pressure_grad_host)
            .with("diff_coeff", diffusion_coeff_host)
            .with("part_flux", particle_flux_host)
            .with("flux_grad", flux_grad_host)
            .with("n_eq_ui", n_eq_ui_host);
}

DFieldSpMomXn DiffGridsFluidSolver::operator()(
        DFieldSpMomXn const neutrals,
        DConstFieldSpXVx const allfdistribu,
        DConstFieldX const efield,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("DiffusiveNeutralSolver");
    RK2<DFieldMemSpMomXn> timestepper(get_idx_range(neutrals));

    // moments computation
    IdxRangeSpX idx_range_kspx(get_idx_range(allfdistribu));
    DFieldMemSpX density_alloc(idx_range_kspx);
    DFieldMemSpX velocity_alloc(idx_range_kspx);
    DFieldMemSpX temperature_alloc(idx_range_kspx);

    DFieldSpX density = get_field(density_alloc);
    DFieldSpX velocity = get_field(velocity_alloc);
    DFieldSpX temperature = get_field(temperature_alloc);

    DConstFieldVx quadrature_coeffs = m_quadrature_coeffs;

    IdxRangeVx const idx_range_vx(get_idx_range<GridVx>(allfdistribu));

    // fluid moments computation
    ddc::parallel_fill(density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_kspx,
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IdxVx const ivx : idx_range_vx) {
                    CoordVx const coordv = ddc::coordinate(ivx);
                    double const val(quadrature_coeffs(ivx) * allfdistribu(ispx, ivx));
                    density(ispx) += val;
                    particle_flux += val * coordv;
                    momentum_flux += val * coordv * coordv;
                }
                velocity(ispx) = particle_flux / density(ispx);
                temperature(ispx)
                        = (momentum_flux - particle_flux * velocity(ispx)) / density(ispx);
            });

    timestepper.update(neutrals, dt, [&](DFieldSpMomXn dn, DConstFieldSpMomXn n) {
        get_derivative(
                dn,
                n,
                get_const_field(density),
                get_const_field(velocity),
                get_const_field(temperature));
    });
    Kokkos::Profiling::popRegion();
    return neutrals;
}

// interpolate to the neutral grid
void DiffGridsFluidSolver::interpolate_on_neutral_grid(
        DFieldSpXn field_on_Xn,
        DConstFieldSpX field_on_X) const
{
    SplineX_GridXnEvaluator interpolator = m_interpolator_from_X_to_Xn;
    ddc::for_each(get_idx_range<Species>(field_on_X), [&](IdxSp const isp) {
        DBSFieldMemX spline_coeff_alloc(get_spline_idx_range(m_spline_builder_on_X));
        DBSFieldX spline_coeff(get_field(spline_coeff_alloc));
        m_spline_builder_on_X(spline_coeff, get_const_field(field_on_X[isp]));

        FieldMemXn<CoordX> coords_eval_alloc(get_idx_range<GridXNeutrals>(field_on_Xn));
        FieldXn<CoordX> coords_eval = get_field(coords_eval_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range<GridXNeutrals>(field_on_Xn),
                KOKKOS_LAMBDA(IdxXn const ixn) { coords_eval(ixn) = ddc::coordinate(ixn); });
        interpolator(field_on_Xn[isp], get_const_field(coords_eval), get_const_field(spline_coeff));
    });
}
