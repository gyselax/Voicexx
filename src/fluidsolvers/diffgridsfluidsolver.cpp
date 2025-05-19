// SPDX-License-Identifier: MIT

#include <string>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "central_fdm_partial_derivatives.hpp"
#include "central_fdm_partial_derivatives_with_boundary_values.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "diffgridsfluidsolver.hpp"
#include "geometry.hpp"
#include "geometry_moments.hpp"
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
        SplineX_GridXnEvaluator const& interpolator_from_X_to_Xn,
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
    , m_interpolator_from_X_to_Xn(interpolator_from_X_to_Xn)
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

void DiffGridsFluidSolver::compute_particle_flux(
        DFieldSpXn const neutrals_particle_flux,
        DConstFieldSpXn const neutrals_density,
        DConstFieldMomSpXn const plasma_moments,
        DConstFieldSpXn const charge_exchange,
        DConstFieldSpXn const ionisation,
        DConstFieldSpXn const recombination) const
{
    DConstFieldSpXn plasma_density = plasma_moments[GeometryMX::density_idx];
    DConstFieldSpXn plasma_velocity = plasma_moments[GeometryMX::velocity_idx];
    DConstFieldSpXn plasma_temperature = plasma_moments[GeometryMX::temperature_idx];

    // compute all the terms needed to compute the particle flux
    IdxRangeSpXn idx_range_neutrals(get_idx_range(neutrals_density));
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

    IdxSp const iion(find_ion(get_idx_range<Species>(plasma_density)));
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
                double density_neutrals = neutrals_density(ifspxn);
                double density_elec = plasma_density(ielec(), ixn);
                double density_ions = plasma_density(iion, ixn);
                double K_cx = charge_exchange(ifspxn);
                double K_i = ionisation(ifspxn);
                double K_r = recombination(ifspxn);
                neutral_pressure(ifspxn) = plasma_temperature(iion, ixn) * density_neutrals;
                density_eff(ifspxn) = (density_neutrals * density_ions * K_cx
                                       + density_ions * density_elec * K_r)
                                      / (density_ions * K_cx + density_elec * K_i);
                density_eff_velocity_ions(ifspxn)
                        = density_eff(ifspxn) * plasma_velocity(iion, ixn);
                pressure_diffusion_coefficient(ifspxn)
                        = mean_free_path / (density_ions * K_cx + density_elec * K_i);
                // density source is not solved here, we only solve transport.
            });

    // the partial derivative creator
    CentralFDMPartialDerivativeCreator<IdxRangeXn, X> const partial_x_creator;

    // compute the gradient of the neutral pressure
    DFieldMemSpXn gradx_pressure_alloc(idx_range_neutrals);
    DFieldSpXn gradx_pressure = get_field(gradx_pressure_alloc);
    ddc::for_each(get_idx_range<Species>(neutrals_density), [&](IdxSp const isp) {
        DFieldXn pressure_sp = neutral_pressure[isp];
        std::unique_ptr<IPartialDerivative<IdxRangeXn, X>> const partial_x_pointer
                = partial_x_creator.create_instance(get_const_field(pressure_sp));
        IPartialDerivative<IdxRangeXn, X> const& partial_x = *partial_x_pointer;
        partial_x(gradx_pressure[isp]);
    });

    // compute the neutral particle flux
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_neutrals,
            KOKKOS_LAMBDA(IdxSpXn const ifspxn) {
                neutrals_particle_flux(ifspxn) = density_eff_velocity_ions(ifspxn) // convection
                                                 - pressure_diffusion_coefficient(ifspxn)
                                                           * gradx_pressure(ifspxn); //diffusion
                neutrals_particle_flux(ifspxn) *= sqrt_mass_ratio;
            });

    // we expose to pdi the coefficients
    auto n_eq_ui_host = ddc::create_mirror_view_and_copy(density_eff_velocity_ions);
    auto pressure_grad_host = ddc::create_mirror_view_and_copy(gradx_pressure);
    auto pressure_host = ddc::create_mirror_view_and_copy(neutral_pressure);
    auto diffusion_coeff_host = ddc::create_mirror_view_and_copy(pressure_diffusion_coefficient);
    ddc::PdiEvent("diff_conv_expose")
            .with("pressure", pressure_host)
            .with("pressure_grad", pressure_grad_host)
            .with("diff_coeff", diffusion_coeff_host)
            .with("n_eq_ui", n_eq_ui_host);
}


void DiffGridsFluidSolver::get_density_derivative(
        DFieldSpXn const derivative_density_neutrals,
        DConstFieldSpXn const density_neutrals,
        DConstFieldMomSpXn const plasma_moments,
        DConstFieldSpXn const charge_exchange,
        DConstFieldSpXn const ionisation,
        DConstFieldSpXn const recombination) const
{
    IdxRangeSpXn idx_range_neutrals(get_idx_range(density_neutrals));
    DFieldMemSpXn particle_flux_alloc(idx_range_neutrals);
    DFieldSpXn particle_flux = get_field(particle_flux_alloc); // \Gamma_N
    compute_particle_flux(
            particle_flux,
            get_const_field(density_neutrals),
            plasma_moments,
            charge_exchange,
            ionisation,
            recombination);

    // compute the divergence of the flux via FDM
    // depending on the boundary condition
    DFieldMemSpXn div_particle_flux_alloc(idx_range_neutrals);
    DFieldSpXn div_particle_flux = get_field(div_particle_flux_alloc); // \partial_x \Gamma_N
    switch (m_flux_BC) {
    case NeutralFluxBoundaryCondition::escaping_neutrals: {
        CentralFDMPartialDerivativeCreator<IdxRangeXn, X> const partial_x_creator;
        ddc::for_each(get_idx_range<Species>(density_neutrals), [&](IdxSp const isp) {
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
        for (IdxSp const isp : get_idx_range<Species>(density_neutrals)) {
            DFieldXn flux_sp = particle_flux[isp];
            std::unique_ptr<IPartialDerivative<IdxRangeXn, X>> const partial_x_pointer
                    = partial_x_bv_creator.create_instance(get_const_field(flux_sp));
            IPartialDerivative<IdxRangeXn, X> const& partial_x = *partial_x_pointer;
            partial_x(div_particle_flux[isp]);
        }
        break;
    }
    case NeutralFluxBoundaryCondition::recycling: {
        IdxRangeXn Xn_range_neutrals(idx_range_neutrals);
        IdxSp const iion(find_ion(get_idx_range<Species>(plasma_moments)));
        IdxXn min(Xn_range_neutrals.front());
        IdxXn max(Xn_range_neutrals.back());
        double const recycling_coeff_proxy = m_recycling_coefficient;
        auto density_proxy = ddc::create_mirror_and_copy(plasma_moments[GeometryMX::density_idx]);
        auto velocity_proxy = ddc::create_mirror_and_copy(plasma_moments[GeometryMX::velocity_idx]);
        double const boundary_condition_left
                = -recycling_coeff_proxy * density_proxy(iion, min) * velocity_proxy(iion, min);
        double const boundary_condition_right
                = -recycling_coeff_proxy * density_proxy(iion, max) * velocity_proxy(iion, max);
        CentralFDMPartialDerivativeWithBValueCreator<IdxRangeXn, X> const
                partial_x_bv_creator(boundary_condition_left, boundary_condition_right);
        for (IdxSp const isp : get_idx_range<Species>(density_neutrals)) {
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

    // Compute the neutral derivative
    // As the equations for the other moments are not solved, we only slice the density
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_neutrals,
            KOKKOS_LAMBDA(IdxSpXn const ifspxn) {
                derivative_density_neutrals(ifspxn) = -div_particle_flux(ifspxn);
            }); // density source is not solved here, we only solve transport.
}

DFieldMomSpXn DiffGridsFluidSolver::operator()(
        DFieldMomSpXn const neutrals_moments,
        DConstFieldSpXVx const allfdistribu,
        DConstFieldX const efield,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("DiffGridsFluidSolver");

    if (!GeometryMX::is_density_and_flux(get_idx_range<GridMom>(neutrals_moments))) {
        throw std::runtime_error(
                "The neutrals fluid moments should contain the density, the fluid velocity, and "
                "nothing more");
    }

    // compute the plasma fluid moments
    IdxRangeSpX idx_range_kinsp_on_X(get_idx_range(allfdistribu));
    IdxRangeMomSpX idx_range_momkspx(GeometryMX::first_three_moments, idx_range_kinsp_on_X);
    DFieldMem<IdxRangeMomSpX> plasma_moments_on_X_alloc(idx_range_momkspx);
    DFieldMomSpX plasma_moments_on_X(plasma_moments_on_X_alloc);
    compute_plasma_moments(plasma_moments_on_X, allfdistribu);

    // interpolate them on the neutrals grid
    IdxRangeMomSpXn idx_range_kinsp_on_Xn(
            GeometryMX::first_three_moments,
            get_idx_range<Species>(allfdistribu),
            get_idx_range<GridXNeutrals>(neutrals_moments));
    DFieldMemMomSpXn plasma_moments_on_Xn_alloc(idx_range_kinsp_on_Xn);
    DFieldMomSpXn plasma_moments_on_Xn(plasma_moments_on_Xn_alloc);
    m_interpolator(
            plasma_moments_on_Xn[GeometryMX::density_idx],
            get_const_field(plasma_moments_on_X[GeometryMX::density_idx]));
    m_interpolator(
            plasma_moments_on_Xn[GeometryMX::velocity_idx],
            get_const_field(plasma_moments_on_X[GeometryMX::velocity_idx]));
    m_interpolator(
            plasma_moments_on_Xn[GeometryMX::temperature_idx],
            get_const_field(plasma_moments_on_X[GeometryMX::temperature_idx]));

    // compute the reaction rates (directly on the neutral grid)
    IdxRangeSpXn idx_range_neutrals(get_idx_range(neutrals_moments));
    DFieldMemSpXn cx_rate_Xn_alloc(idx_range_neutrals);
    DFieldMemSpXn i_rate_Xn_alloc(idx_range_neutrals);
    DFieldMemSpXn r_rate_Xn_alloc(idx_range_neutrals);
    DFieldSpXn charge_exchange_rate_on_Xn = get_field(cx_rate_Xn_alloc);
    DFieldSpXn ionisation_rate_on_Xn = get_field(i_rate_Xn_alloc);
    DFieldSpXn recombination_rate_on_Xn = get_field(r_rate_Xn_alloc);
    compute_reaction_rates(
            charge_exchange_rate_on_Xn,
            ionisation_rate_on_Xn,
            recombination_rate_on_Xn,
            get_const_field(neutrals_moments),
            get_const_field(plasma_moments_on_X));

    // we store in neutrals_moments the particle flux, which is neither solved using
    // an equation nor used here. We store it to give it to the plasma-neutrals
    // coupling operator that needs it afterwards.
    DFieldSpXn neutrals_density = neutrals_moments[GeometryMX::density_idx];
    DFieldSpXn neutrals_particle_flux = neutrals_moments[GeometryMX::velocity_idx];
    compute_particle_flux(
            neutrals_particle_flux,
            get_const_field(neutrals_density),
            get_const_field(plasma_moments_on_Xn),
            get_const_field(charge_exchange_rate_on_Xn),
            get_const_field(ionisation_rate_on_Xn),
            get_const_field(recombination_rate_on_Xn));

    // we only update the density because it is the only moment we solve
    RK2<DFieldMemSpXn> timestepper(get_idx_range(neutrals_density));
    timestepper.update(neutrals_density, dt, [&](DFieldSpXn dn, DConstFieldSpXn n) {
        get_density_derivative(
                dn,
                n,
                get_const_field(plasma_moments_on_Xn),
                get_const_field(charge_exchange_rate_on_Xn),
                get_const_field(ionisation_rate_on_Xn),
                get_const_field(recombination_rate_on_Xn));
    });

    Kokkos::Profiling::popRegion();
    return neutrals_moments;
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

void DiffGridsFluidSolver::compute_plasma_moments(
        DFieldMomSpX const plasma_moments,
        DConstFieldSpXVx const allfdistribu) const
{
    IdxRangeSpX idx_range_kspx(get_idx_range(allfdistribu));
    DFieldSpX density = plasma_moments[GeometryMX::density_idx];
    DFieldSpX velocity = plasma_moments[GeometryMX::velocity_idx];
    DFieldSpX temperature = plasma_moments[GeometryMX::temperature_idx];

    DConstFieldVx quadrature_coeffs = m_quadrature_coeffs;

    IdxRangeVx const idx_range_vx(get_idx_range<GridVx>(allfdistribu));

    // fluid moments computation
    ddc::parallel_fill(plasma_moments, 0.);
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
}

void DiffGridsFluidSolver::compute_reaction_rates(
        DFieldSpXn const charge_exchange_on_Xn,
        DFieldSpXn const ionisation_on_Xn,
        DFieldSpXn const recombination_on_Xn,
        DConstFieldMomSpXn const neutrals_moments,
        DConstFieldMomSpX const plasma_moments_on_X) const
{
    IdxRangeSpX idx_range_rates_on_X(
            get_idx_range<Species>(neutrals_moments),
            get_idx_range<GridX>(plasma_moments_on_X));

    // building reaction rates
    // pay attention, these are normalised to Kcx0=10^-14
    DFieldMemSpX charge_exchange_rate_alloc(idx_range_rates_on_X);
    DFieldMemSpX ionisation_rate_alloc(idx_range_rates_on_X);
    DFieldMemSpX recombination_rate_alloc(idx_range_rates_on_X);

    DFieldSpX charge_exchange_rate = get_field(charge_exchange_rate_alloc);
    DFieldSpX ionisation_rate = get_field(ionisation_rate_alloc);
    DFieldSpX recombination_rate = get_field(recombination_rate_alloc);

    DConstFieldSpX plasma_density = plasma_moments_on_X[GeometryMX::density_idx];
    DConstFieldSpX plasma_temperature = plasma_moments_on_X[GeometryMX::temperature_idx];

    m_charge_exchange(charge_exchange_rate, plasma_density, plasma_temperature);
    m_ionisation(ionisation_rate, plasma_density, plasma_temperature);
    m_recombination(recombination_rate, plasma_density, plasma_temperature);

    // expose to pdi the reaction rate coefficients
    auto cx_host = ddc::create_mirror_view_and_copy(charge_exchange_rate);
    auto i_host = ddc::create_mirror_view_and_copy(ionisation_rate);
    auto r_host = ddc::create_mirror_view_and_copy(recombination_rate);
    ddc::PdiEvent("reaction_rate_expose")
            .with("charge_exchange_rate", cx_host)
            .with("ionisation_rate", i_host)
            .with("recombination_rate", r_host);

    m_interpolator(charge_exchange_on_Xn, get_const_field(charge_exchange_rate));
    m_interpolator(ionisation_on_Xn, get_const_field(ionisation_rate));
    m_interpolator(recombination_on_Xn, get_const_field(recombination_rate));
}
