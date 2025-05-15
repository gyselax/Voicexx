// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "mask_tanh.hpp"
#include "rk2.hpp"
#include "samegridfluidsolver.hpp"

SameGridFluidSolver::SameGridFluidSolver(
        IReactionRate const& charge_exchange,
        IReactionRate const& ionisation,
        IReactionRate const& recombination,
        double const normalisation_coeff,
        SplineXBuilder const& spline_x_builder,
        SplineXEvaluator const& spline_x_evaluator,
        DConstFieldVx const& quadrature_coeffs,
        double const neutrals_wall_extent,
        double const neutrals_wall_stiffness,
        double const neutrals_wall_amplitude,
        IdxRangeX const& gridx)

    : m_charge_exchange(charge_exchange)
    , m_ionisation(ionisation)
    , m_recombination(recombination)
    , m_normalisation_coeff(normalisation_coeff)
    , m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_quadrature_coeffs(quadrature_coeffs)
    , m_mask_amplitude(neutrals_wall_amplitude)
    , m_mask(gridx)
{
    // The amplitude cannot be negative
    if (m_mask_amplitude < 0) {
        throw std::runtime_error("The amplitude should be positive");
    }

    host_t<DFieldMemX> mask_host(gridx);
    mask_host = mask_tanh(
            gridx,
            neutrals_wall_extent,
            neutrals_wall_stiffness,
            MaskType::Inverted,
            false);
    ddc::parallel_deepcopy(get_field(m_mask), mask_host);
    ddc::expose_to_pdi("krook_neutrals_amplitude", neutrals_wall_amplitude);
    ddc::expose_to_pdi("krook_neutrals_mask", mask_host);
}

IdxSp SameGridFluidSolver::find_ion(IdxRangeSp const idx_range_kinsp) const
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

void SameGridFluidSolver::get_derivative(
        DFieldSpMomX dn,
        DConstFieldSpMomX neutrals,
        DConstFieldSpX density,
        DConstFieldSpX velocity,
        DConstFieldSpX temperature) const
{
    IdxRangeSpX idx_range_fluidspx(get_idx_range<Species, GridX>(neutrals));

    // building reaction rates
    DFieldMemSpX charge_exchange_rate_alloc(idx_range_fluidspx);
    DFieldMemSpX ionisation_rate_alloc(idx_range_fluidspx);
    DFieldMemSpX recombination_rate_alloc(idx_range_fluidspx);

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

    // compute diffusive model equation terms
    DFieldMemSpX density_equilibrium_velocity_alloc(idx_range_fluidspx);
    DFieldMemSpX diffusion_temperature_alloc(idx_range_fluidspx);
    DFieldSpX density_equilibrium_velocity = get_field(density_equilibrium_velocity_alloc);
    DFieldSpX diffusion_temperature = get_field(diffusion_temperature_alloc);

    IdxSp const iion(find_ion(get_idx_range<Species>(density)));
    IdxMom const ineutral_density(0);

    double const normalisation_coeff_alpha0(m_normalisation_coeff);
    double const mass_ratio(mass(ielec()) / mass(iion));
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_fluidspx,
            KOKKOS_LAMBDA(IdxSpX const ifspx) {
                IdxSp const isp(ddc::select<Species>(ifspx));
                IdxX const ix(ddc::select<GridX>(ifspx));

                double const denom = density(iion, ix) * charge_exchange_rate(ifspx)
                                     + density(ielec(), ix) * ionisation_rate(ifspx);

                density_equilibrium_velocity(ifspx)
                        = (density(ielec(), ix) * density(iion, ix) * recombination_rate(ifspx)
                           + neutrals(ifspx, ineutral_density) * density(iion, ix)
                                     * charge_exchange_rate(ifspx))
                          * velocity(iion, ix) * Kokkos::sqrt(mass_ratio) / denom;

                diffusion_temperature(ifspx)
                        = normalisation_coeff_alpha0 * temperature(iion, ix) / (mass(isp) * denom);
                // density source is not solved here, we only solve transport.
            });

    // compute coordinates at which spatial derivatives are evaluated
    FieldMemX<CoordX> coords_eval_alloc(get_idx_range<GridX>(neutrals));
    FieldX<CoordX> coords_eval = get_field(coords_eval_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range<GridX>(neutrals),
            KOKKOS_LAMBDA(IdxX const ix) { coords_eval(ix) = ddc::coordinate(ix); });

    // create chunks to store spatial derivatives
    DFieldMemSpX gradx_density_equilibrium_velocity_alloc(idx_range_fluidspx);
    DFieldMemSpX gradx_diffusion_temperature_alloc(idx_range_fluidspx);
    DFieldMemSpX gradx_neutrals_density_alloc(idx_range_fluidspx);
    DFieldMemSpX laplx_neutrals_density_alloc(idx_range_fluidspx);

    DFieldSpX gradx_density_equilibrium_velocity
            = get_field(gradx_density_equilibrium_velocity_alloc);
    DFieldSpX gradx_diffusion_temperature = get_field(gradx_diffusion_temperature_alloc);
    DFieldSpX gradx_neutrals_density = get_field(gradx_neutrals_density_alloc);
    DFieldSpX laplx_neutrals_density = get_field(laplx_neutrals_density_alloc);

    ddc::for_each(get_idx_range<Species>(neutrals), [&](IdxSp const isp) {
        // compute spline coefficients
        DBSFieldMemX density_equilibrium_velocity_spline_x_coeff(
                get_spline_idx_range(m_spline_x_builder));

        DBSFieldMemX diffusion_temperature_spline_x_coeff(get_spline_idx_range(m_spline_x_builder));

        DBSFieldMemX neutrals_density_spline_x_coeff(get_spline_idx_range(m_spline_x_builder));

        m_spline_x_builder(
                get_field(density_equilibrium_velocity_spline_x_coeff),
                get_const_field(density_equilibrium_velocity[isp]));

        m_spline_x_builder(
                get_field(diffusion_temperature_spline_x_coeff),
                get_const_field(diffusion_temperature[isp]));

        m_spline_x_builder(
                get_field(neutrals_density_spline_x_coeff),
                get_const_field(neutrals[IdxSpMom(isp, ineutral_density)]));

        // compute gradients
        m_spline_x_evaluator
                .deriv(get_field(gradx_density_equilibrium_velocity[isp]),
                       get_const_field(coords_eval),
                       get_const_field(density_equilibrium_velocity_spline_x_coeff));

        m_spline_x_evaluator
                .deriv(gradx_diffusion_temperature[isp],
                       get_const_field(coords_eval),
                       get_const_field(diffusion_temperature_spline_x_coeff));

        m_spline_x_evaluator
                .deriv(gradx_neutrals_density[isp],
                       get_const_field(coords_eval),
                       get_const_field(neutrals_density_spline_x_coeff));

        // compute laplacian
        DBSFieldMemX gradx_neutrals_density_spline_x_coeff(
                get_spline_idx_range(m_spline_x_builder));

        m_spline_x_builder(
                get_field(gradx_neutrals_density_spline_x_coeff),
                get_const_field(gradx_neutrals_density[isp]));

        m_spline_x_evaluator
                .deriv(laplx_neutrals_density[isp],
                       get_const_field(coords_eval),
                       get_const_field(gradx_neutrals_density_spline_x_coeff));
    });

    // get the neutral mask
    DConstFieldX mask(get_field(m_mask));
    double amplitude = m_mask_amplitude;

    DFieldMemSpX diff_term_alloc(idx_range_fluidspx);
    DFieldMemSpX conv_term_alloc(idx_range_fluidspx);
    DFieldSpX diff_term = get_field(diff_term_alloc);
    DFieldSpX conv_term = get_field(conv_term_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_fluidspx,
            KOKKOS_LAMBDA(IdxSpX const ifspx) {
                IdxX ifx(ifspx);
                diff_term(ifspx) = -gradx_density_equilibrium_velocity(ifspx);
                conv_term(ifspx)
                        = gradx_diffusion_temperature(ifspx) * gradx_neutrals_density(ifspx)
                          + diffusion_temperature(ifspx) * laplx_neutrals_density(ifspx);
                dn(ifspx, ineutral_density)
                        // here the masks stops the diffusion and convection in the wall
                        = (diff_term(ifspx) + conv_term(ifspx)) * (1. - mask(ifx))
                          // we relax the density to 0 in the walls
                          - amplitude * (mask(ifx)) * neutrals(ifspx, ineutral_density);
            }); // density source is not solved here, we only solve transport.

    // we expose to pdi the coefficients
    auto diff_term_host = ddc::create_mirror_view_and_copy(diff_term);
    auto conv_term_host = ddc::create_mirror_view_and_copy(conv_term);
    ddc::PdiEvent("diff_conv_expose")
            .with("diffusion_term", diff_term_host)
            .with("convection_term", conv_term_host);
}

DFieldSpMomX SameGridFluidSolver::operator()(
        DFieldSpMomX const neutrals,
        DConstFieldSpXVx const allfdistribu,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("SameGridFluidSolver");
    RK2<DFieldMemSpMomX> timestepper(get_idx_range(neutrals));

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

    timestepper.update(neutrals, dt, [&](DFieldSpMomX dn, DConstFieldSpMomX n) {
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
