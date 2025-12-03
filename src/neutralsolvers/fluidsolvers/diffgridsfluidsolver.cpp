// SPDX-License-Identifier: MIT

#include <string>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "central_fdm_partial_derivatives.hpp"
#include "central_fdm_partial_derivatives_with_boundary_values.hpp"
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
        GridNeutralInterpolator const& interpolator,
        DConstFieldVx const& quadrature_coeffs,
        NeutralFluxBoundaryCondition flux_BC,
        double recycling_coeff)

    : m_charge_exchange(charge_exchange)
    , m_ionisation(ionisation)
    , m_recombination(recombination)
    , m_mean_free_path(mean_free_path)
    , m_spline_builder_on_Xn(spline_builder_on_Xn)
    , m_spline_evaluator_on_Xn(spline_evaluator_on_Xn)
    , m_spline_builder_on_X(spline_builder_on_X)
    , m_interpolator(interpolator)
    , m_quadrature_coeffs(quadrature_coeffs)
    , m_flux_BC(flux_BC)
    , m_recycling_coefficient(recycling_coeff)
{
}

IdxSp DiffGridsFluidSolver::find_ion(IdxRangeSp const idx_range_kinsp) const
{
    bool ion_found = false;
    IdxSp i_ion(0);
    for (IdxSp const isp : idx_range_kinsp) {
        if (charge(isp) > 0.) {
            ion_found = true;
            i_ion = isp;
        }
    }
    if (!ion_found) {
        throw std::runtime_error("ion not found");
    }
    assert(idx_range_kinsp.size() == 2);

    return i_ion;
}

void DiffGridsFluidSolver::get_derivative(
        DFieldSpMomXn dn,
        DConstFieldSpMomXn neutrals,
        DConstFieldSpX density,
        DConstFieldSpX velocity,
        DConstFieldSpX temperature) const
{
    IdxRangeSpXn idx_range_neutrals(get_idx_range(neutrals));
    IdxRangeSpX
            idx_range_rates_on_X(get_idx_range<Species>(neutrals), get_idx_range<GridX>(density));

    // building reaction rates
    // pay attention, these are normalised to Kcx0=10^-14
    DFieldMemSpX charge_exchange_rate_alloc(idx_range_rates_on_X);
    DFieldMemSpX ionisation_rate_alloc(idx_range_rates_on_X);
    DFieldMemSpX recombination_rate_alloc(idx_range_rates_on_X);

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
    IdxRangeSpXn idx_range_kinsp_onXn(
            get_idx_range<Species>(density),
            get_idx_range<GridXNeutrals>(neutrals));
    DFieldMemSpXn density_on_Xn_alloc(idx_range_kinsp_onXn);
    DFieldSpXn density_on_Xn = get_field(density_on_Xn_alloc);
    DFieldMemSpXn velocity_on_Xn_alloc(idx_range_kinsp_onXn);
    DFieldSpXn velocity_on_Xn = get_field(velocity_on_Xn_alloc);
    DFieldMemSpXn temperature_on_Xn_alloc(idx_range_kinsp_onXn);
    DFieldSpXn temperature_on_Xn = get_field(temperature_on_Xn_alloc);

    DFieldMemSpXn cx_rate_Xn_alloc(idx_range_neutrals);
    DFieldMemSpXn i_rate_Xn_alloc(idx_range_neutrals);
    DFieldMemSpXn r_rate_Xn_alloc(idx_range_neutrals);
    DFieldSpXn charge_exchange_rate_on_Xn = get_field(cx_rate_Xn_alloc);
    DFieldSpXn ionisation_rate_on_Xn = get_field(i_rate_Xn_alloc);
    DFieldSpXn recombination_rate_on_Xn = get_field(r_rate_Xn_alloc);

    // do the interpolation on the neutral grid
    m_interpolator(density_on_Xn, density);
    m_interpolator(velocity_on_Xn, velocity);
    m_interpolator(temperature_on_Xn, temperature);

    m_interpolator(charge_exchange_rate_on_Xn, get_const_field(charge_exchange_rate));
    m_interpolator(ionisation_rate_on_Xn, get_const_field(ionisation_rate));
    m_interpolator(recombination_rate_on_Xn, get_const_field(recombination_rate));

    // compute diffusive model equation terms
    DFieldMemSpXn density_eff_alloc(idx_range_neutrals);
    DFieldMemSpXn density_eff_velocity_ions_alloc(idx_range_neutrals);
    DFieldMemSpXn neutral_pressure_alloc(idx_range_neutrals);
    DFieldMemSpXn pressure_diffusion_coefficient_alloc(idx_range_neutrals);
    DFieldSpXn density_eff = get_field(density_eff_alloc); // n_eff
    DFieldSpXn density_eff_velocity_ions
            = get_field(density_eff_velocity_ions_alloc); // n_eff * u_i
    DFieldSpXn neutral_pressure
            = get_field(neutral_pressure_alloc); // p_N=T_i*n_N because (T_i=T_N)
    DFieldSpXn pressure_diffusion_coefficient
            = get_field(pressure_diffusion_coefficient_alloc); // \hat D

    IdxSp const iion(find_ion(get_idx_range<Species>(density)));
    IdxMom const ineutral_density(0);
    double const sqrt_mass_ratio(Kokkos::sqrt(mass(ielec()) / mass(iion)));
    double const mean_free_path = m_mean_free_path;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_neutrals,
            KOKKOS_LAMBDA(IdxSpXn const ifspxn) {
                IdxSp const isp(ifspxn);
                IdxXn const ixn(ifspxn);

                // in this for loop we construct the terms n_eff*u_i, p_N and \hat D
                // we will take their derivatives afterwards
                double density_neutrals = neutrals(ifspxn, ineutral_density);
                double density_elec = density_on_Xn(ielec(), ixn);
                double density_ions = density_on_Xn(iion, ixn);
                double K_cx = charge_exchange_rate_on_Xn(ifspxn);
                double K_i = ionisation_rate_on_Xn(ifspxn);
                double K_r = recombination_rate_on_Xn(ifspxn);
                neutral_pressure(ifspxn) = temperature_on_Xn(iion, ixn) * density_neutrals;
                density_eff(ifspxn) = (density_neutrals * density_ions * K_cx
                                       + density_ions * density_elec * K_r)
                                      / (density_ions * K_cx + density_elec * K_i);
                density_eff_velocity_ions(ifspxn) = density_eff(ifspxn) * velocity_on_Xn(iion, ixn);
                pressure_diffusion_coefficient(ifspxn)
                        = mean_free_path / (density_ions * K_cx + density_elec * K_i);
                // density source is not solved here, we only solve transport.
            });

    // the partial derivative creator
    CentralFDMPartialDerivativeCreator<IdxRangeXn, X> const partial_x_creator;

    // compute the gradient of the neutral pressure
    IdxRangeXn Xn_range_neutrals(idx_range_neutrals);
    DFieldMemSpXn gradx_pressure_alloc(idx_range_neutrals);
    DFieldSpXn gradx_pressure = get_field(gradx_pressure_alloc);
    ddc::for_each(get_idx_range<Species>(neutrals), [&](IdxSp const isp) {
        DFieldXn pressure_sp = neutral_pressure[isp];
        std::unique_ptr<IPartialDerivative<IdxRangeXn, X>> const partial_x_pointer
                = partial_x_creator.create_instance(get_const_field(pressure_sp));
        IPartialDerivative<IdxRangeXn, X> const& partial_x = *partial_x_pointer;
        partial_x(gradx_pressure[isp]);
    });

    // compute the neutral particle flux
    DFieldMemSpXn particle_flux_alloc(idx_range_neutrals);
    DFieldSpXn particle_flux = get_field(particle_flux_alloc); // \Gamma_N
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_neutrals,
            KOKKOS_LAMBDA(IdxSpXn const ifspxn) {
                particle_flux(ifspxn) = density_eff_velocity_ions(ifspxn) // convection
                                        - pressure_diffusion_coefficient(ifspxn)
                                                  * gradx_pressure(ifspxn); //diffusion
            });

    // compute the divergence of the flux via FDM
    // depending on the boundary condition
    DFieldMemSpXn div_particle_flux_alloc(idx_range_neutrals);
    DFieldSpXn div_particle_flux = get_field(div_particle_flux_alloc); // \partial_x \Gamma_N
    switch (m_flux_BC) {
    case NeutralFluxBoundaryCondition::escaping_neutrals: {
        ddc::for_each(get_idx_range<Species>(neutrals), [&](IdxSp const isp) {
            DFieldXn flux_sp = particle_flux[isp];
            std::unique_ptr<IPartialDerivative<IdxRangeXn, X>> const partial_x_pointer
                    = partial_x_creator.create_instance(get_const_field(flux_sp));
            IPartialDerivative<IdxRangeXn, X> const& partial_x = *partial_x_pointer;
            partial_x(div_particle_flux[isp]);
        });
        break;
    }
    case NeutralFluxBoundaryCondition::zero_flux: {
        double const boundary_condition_flux = 0;
        CentralFDMPartialDerivativeWithBValueCreator<IdxRangeXn, X> const
                partial_x_bv_creator(boundary_condition_flux, boundary_condition_flux);
        for (IdxSp const isp : get_idx_range<Species>(neutrals)) {
            DFieldXn flux_sp = particle_flux[isp];
            std::unique_ptr<IPartialDerivative<IdxRangeXn, X>> const partial_x_pointer
                    = partial_x_bv_creator.create_instance(get_const_field(flux_sp));
            IPartialDerivative<IdxRangeXn, X> const& partial_x = *partial_x_pointer;
            partial_x(div_particle_flux[isp]);
        }
        break;
    }
    case NeutralFluxBoundaryCondition::recycling: {
        IdxXn min(Xn_range_neutrals.front());
        IdxXn max(Xn_range_neutrals.back());
        double const recycling_coeff_proxy = m_recycling_coefficient;
        auto density_proxy = ddc::create_mirror_and_copy(density_on_Xn);
        auto velocity_proxy = ddc::create_mirror_and_copy(velocity_on_Xn);
        double const boundary_condition_left
                = -recycling_coeff_proxy * density_proxy(iion, min) * velocity_proxy(iion, min);
        double const boundary_condition_right
                = -recycling_coeff_proxy * density_proxy(iion, max) * velocity_proxy(iion, max);
        CentralFDMPartialDerivativeWithBValueCreator<IdxRangeXn, X> const
                partial_x_bv_creator(boundary_condition_left, boundary_condition_right);
        for (IdxSp const isp : get_idx_range<Species>(neutrals)) {
            DFieldXn flux_sp = particle_flux[isp];
            std::unique_ptr<IPartialDerivative<IdxRangeXn, X>> const partial_x_pointer
                    = partial_x_bv_creator.create_instance(get_const_field(flux_sp));
            IPartialDerivative<IdxRangeXn, X> const& partial_x = *partial_x_pointer;
            partial_x(div_particle_flux[isp]);
        }
        break;
    }
    default: {
        throw std::runtime_error("This type of boundary condition is not implemented");
        break;
    }
    }

    // compute the neutral derivative
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_neutrals,
            KOKKOS_LAMBDA(IdxSpXn const ifspxn) {
                dn(ifspxn, ineutral_density) = -sqrt_mass_ratio * div_particle_flux(ifspxn);
            }); // density source is not solved here, we only solve transport.

    // we expose to pdi the coefficients
    auto n_eq_ui_host = ddc::create_mirror_view_and_copy(density_eff_velocity_ions);
    auto pressure_grad_host = ddc::create_mirror_view_and_copy(gradx_pressure);
    auto pressure_host = ddc::create_mirror_view_and_copy(neutral_pressure);
    auto diffusion_coeff_host = ddc::create_mirror_view_and_copy(pressure_diffusion_coefficient);
    auto particle_flux_host = ddc::create_mirror_view_and_copy(particle_flux);
    auto flux_div_host = ddc::create_mirror_view_and_copy(div_particle_flux);
    ddc::PdiEvent("diff_conv_expose")
            .with("pressure", pressure_host)
            .with("pressure_grad", pressure_grad_host)
            .with("diff_coeff", diffusion_coeff_host)
            .with("part_flux", particle_flux_host)
            .with("flux_grad", flux_div_host)
            .with("n_eq_ui", n_eq_ui_host);
}

DFieldSpMomXn DiffGridsFluidSolver::operator()(
        DFieldSpMomXn const neutrals,
        DConstFieldSpXVx const allfdistribu,
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

DiffGridsFluidSolver::NeutralFluxBoundaryCondition DiffGridsFluidSolver::
        neutral_flux_boundary_condition(std::string const& bc_flux_input)
{
    DiffGridsFluidSolver::NeutralFluxBoundaryCondition bc_flux;
    if (bc_flux_input == "recycling") {
        bc_flux = DiffGridsFluidSolver::NeutralFluxBoundaryCondition::recycling;
    } else if (bc_flux_input == "escaping neutrals") {
        bc_flux = DiffGridsFluidSolver::NeutralFluxBoundaryCondition::escaping_neutrals;
    } else if (bc_flux_input == "zero flux") {
        bc_flux = DiffGridsFluidSolver::NeutralFluxBoundaryCondition::zero_flux;
    } else {
        bc_flux = DiffGridsFluidSolver::NeutralFluxBoundaryCondition::escaping_neutrals;
        printf("WARNING: No valid boundary condition for the neutral flux, escaping neutrals "
               "has been taken.\n");
    }
    return bc_flux;
}
