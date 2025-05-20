// SPDX-License-Identifier: MIT

#include <stdexcept>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "densitycoupling.hpp"
#include "geometry.hpp"
#include "geometry_moments.hpp"
#include "geometry_neutrals.hpp"
#include "gridneutral_interpolator.hpp"
#include "ireactionrate.hpp"
#include "rk2.hpp"
#include "species_info.hpp"

DensityCoupling::DensityCoupling(
        double const density_coupling_coeff,
        double const momentum_coupling_coeff,
        double const energy_coupling_coeff,
        IReactionRate const& ionisation,
        IReactionRate const& recombination,
        SplineXBuilder const& spline_builder_on_X,
        GridNeutralInterpolator const& interpolator_between_X_and_Xn,
        double const mean_free_path,
        DConstFieldVx const& quadrature_coeffs)
    : m_density_coupling_coeff(density_coupling_coeff)
    , m_momentum_coupling_coeff(momentum_coupling_coeff)
    , m_energy_coupling_coeff(energy_coupling_coeff)
    , m_ionisation(ionisation)
    , m_recombination(recombination)
    , m_mean_free_path(mean_free_path)
    , m_quadrature_coeffs(quadrature_coeffs)
    , m_spline_builder_on_X(spline_builder_on_X)
    , m_interpolator_between_X_and_Xn(interpolator_between_X_and_Xn)
{
    ddc::expose_to_pdi(
            "kinetic_fluid_coupling_source_density_coupling_coeff",
            m_density_coupling_coeff);
    ddc::expose_to_pdi(
            "kinetic_fluid_coupling_source_momentum_coupling_coeff",
            m_momentum_coupling_coeff);
    ddc::expose_to_pdi(
            "kinetic_fluid_coupling_source_energy_coupling_coeff",
            m_energy_coupling_coeff);
}

IdxSp DensityCoupling::find_ion(IdxRangeSp const dom_kinsp) const
{
    bool ion_found = false;
    IdxSp iion;
    for (IdxSp const isp : dom_kinsp) {
        if (charge(isp) > 0) {
            ion_found = true;
            iion = isp;
        }
    }
    if (!ion_found) {
        throw std::runtime_error("ion not found");
    }
    assert(dom_kinsp.size() == 2);
    return iion;
}

void DensityCoupling::get_particle_source_term(
        DFieldSpXn density_source_neutral,
        DConstFieldSpXn kinsp_density,
        DConstFieldSpXn density_neutrals,
        DConstFieldSpXn ionisation,
        DConstFieldSpXn recombination) const
{
    IdxSp const iion(find_ion(get_idx_range<Species>(kinsp_density)));
    double const sqrt_mass_ratio(Kokkos::sqrt(mass(ielec()) / mass(iion)));
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(density_source_neutral),
            KOKKOS_LAMBDA(IdxSpXn const ispxn) {
                IdxXn ixn(ispxn);
                density_source_neutral(ispxn)
                        = -density_neutrals(ispxn) * kinsp_density(ielec(), ixn) * ionisation(ispxn)
                          + kinsp_density(iion, ixn) * kinsp_density(ielec(), ixn)
                                    * recombination(ispxn);
                density_source_neutral(ispxn) *= sqrt_mass_ratio;
            });
}

void DensityCoupling::get_plasma_source_term(
        DFieldSpXVx plasma_source_term,
        DConstFieldSpX kinsp_temperature,
        DConstFieldSpXn neutral_particle_source_on_Xn) const
{
    // interpolate the particle source on GridX
    IdxRangeSp neutrals_species(get_idx_range<Species>(neutral_particle_source_on_Xn));
    if (neutrals_species.size() != 1) {
        throw std::runtime_error(
                "For the moments the coupling operator only works for one neutrals species");
    }
    IdxRangeSpX
            idx_range_particle_source(neutrals_species, get_idx_range<GridX>(plasma_source_term));
    DFieldMemSpX particle_source_alloc(idx_range_particle_source);
    DFieldSpX particle_source = get_field(particle_source_alloc);
    m_interpolator_between_X_and_Xn(particle_source, neutral_particle_source_on_Xn);
    DConstFieldX plasma_particle_source
            = get_const_field(particle_source[neutrals_species.front()]);

    /*double density_coupling_coeff_proxy = m_density_coupling_coeff;*/
    /*double momentum_coupling_coeff_proxy = m_momentum_coupling_coeff;*/
    /*double energy_coupling_coeff_proxy = m_energy_coupling_coeff;*/
    /*double mean_free_path_proxy = m_mean_free_path;*/

    // ion species
    IdxSp const iion(find_ion(get_idx_range<Species>(plasma_source_term)));
    double const mass_ion = mass(iion);
    DField<IdxRangeXVx> ions_source_term = plasma_source_term[iion];
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(ions_source_term),
            KOKKOS_LAMBDA(IdxXVx const ixvx) {
                ions_source_term(ixvx) = 0;
                IdxX const ix(ixvx);
                IdxVx const ivx(ixvx);
                CoordVx const coordvx = ddc::coordinate(ivx);
                double const density_source = plasma_particle_source(ix); // S_{n,N}
                double const momentum_source = 0; // S_{m,i}
                double const velocity_source = momentum_source / density_source; // U_{N,i}
                double const energy_source = 0; // S_{E,i}
                double const temperature_source = 2 * (energy_source / density_source)
                                                  - mass(iion) * velocity_source; // T_{N,i}
                /*double const neutral_temperature = kinsp_temperature(iion, ix);*/
                /*double const coordvx_sq = coordvx * coordvx;*/
                /*double const density_source*/
                /*        = density_coupling_coeff_proxy*/
                /*          * (1.5 - coordvx_sq / (2 * neutral_temperature))*/
                /*          * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));*/
                /*double const momentum_source*/
                /*        = momentum_coupling_coeff_proxy * Kokkos::sqrt(2 / neutral_temperature)*/
                /*          * coordvx * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));*/
                /*double const energy_source = 2 * energy_coupling_coeff_proxy*/
                /*                             * (-1 + coordvx_sq / neutral_temperature)*/
                /*                             * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));*/
                ions_source_term(ixvx)
                        = -(density_source
                            / (Kokkos::sqrt(2 * M_PI * temperature_source / mass_ion)))
                          * Kokkos::exp(
                                  -(mass_ion * Kokkos::pow(coordvx - velocity_source, 2))
                                  / (2 * temperature_source));
            });

    // electron species
    DField<IdxRangeXVx> electrons_source_term = plasma_source_term[ielec()];
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(ions_source_term),
            KOKKOS_LAMBDA(IdxXVx const ixvx) { electrons_source_term(ixvx) = 0; });
}

void DensityCoupling::get_derivative_neutrals(
        DFieldSpXn dn,
        DConstFieldSpXn particle_source_neutrals) const
{
    IdxRangeSpXn range_neutrals_spxn(get_idx_range(particle_source_neutrals));
    double mean_free_path_proxy = m_mean_free_path;
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            range_neutrals_spxn,
            KOKKOS_LAMBDA(IdxSpXn const ifspx) {
                dn(ifspx) = particle_source_neutrals(ifspx) / mean_free_path_proxy;
            });
}

void DensityCoupling::get_derivative_allfdistribu(
        DFieldSpXVx df,
        DConstFieldSpXVx plasma_source_term) const
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(plasma_source_term),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) { df(ispxvx) = plasma_source_term(ispxvx); });
}

void DensityCoupling::compute_plasma_moments(
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

void DensityCoupling::compute_reaction_rates(
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
    DFieldMemSpX ionisation_rate_alloc(idx_range_rates_on_X);
    DFieldMemSpX recombination_rate_alloc(idx_range_rates_on_X);

    DFieldSpX ionisation_rate = get_field(ionisation_rate_alloc);
    DFieldSpX recombination_rate = get_field(recombination_rate_alloc);

    DConstFieldSpX plasma_density = plasma_moments_on_X[GeometryMX::density_idx];
    DConstFieldSpX plasma_temperature = plasma_moments_on_X[GeometryMX::temperature_idx];

    m_ionisation(ionisation_rate, plasma_density, plasma_temperature);
    m_recombination(recombination_rate, plasma_density, plasma_temperature);

    m_interpolator_between_X_and_Xn(ionisation_on_Xn, get_const_field(ionisation_rate));
    m_interpolator_between_X_and_Xn(recombination_on_Xn, get_const_field(recombination_rate));
}

void DensityCoupling::operator()(
        DFieldSpXVx const allfdistribu,
        DFieldMomSpXn neutrals_moments,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("KineticFluidCouplingSource");

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
    compute_plasma_moments(plasma_moments_on_X, get_const_field(allfdistribu));

    // interpolate them on the neutrals grid
    IdxRangeMomSpXn idx_range_kinsp_on_Xn(
            GeometryMX::first_three_moments,
            get_idx_range<Species>(allfdistribu),
            get_idx_range<GridXNeutrals>(neutrals_moments));
    DFieldMemMomSpXn plasma_moments_on_Xn_alloc(idx_range_kinsp_on_Xn);
    DFieldMomSpXn plasma_moments_on_Xn(plasma_moments_on_Xn_alloc);
    m_interpolator_between_X_and_Xn(
            plasma_moments_on_Xn[GeometryMX::density_idx],
            get_const_field(plasma_moments_on_X[GeometryMX::density_idx]));
    m_interpolator_between_X_and_Xn(
            plasma_moments_on_Xn[GeometryMX::velocity_idx],
            get_const_field(plasma_moments_on_X[GeometryMX::velocity_idx]));
    m_interpolator_between_X_and_Xn(
            plasma_moments_on_Xn[GeometryMX::temperature_idx],
            get_const_field(plasma_moments_on_X[GeometryMX::temperature_idx]));

    // compute the reaction rates (directly on the neutral grid)
    IdxRangeSpXn idx_range_neutrals(get_idx_range(neutrals_moments));
    DFieldMemSpXn i_rate_Xn_alloc(idx_range_neutrals);
    DFieldMemSpXn r_rate_Xn_alloc(idx_range_neutrals);
    DFieldSpXn ionisation_rate_on_Xn = get_field(i_rate_Xn_alloc);
    DFieldSpXn recombination_rate_on_Xn = get_field(r_rate_Xn_alloc);
    compute_reaction_rates(
            ionisation_rate_on_Xn,
            recombination_rate_on_Xn,
            get_const_field(neutrals_moments),
            get_const_field(plasma_moments_on_X));

    // particle source term computation, on Xn
    DFieldMemSpXn particle_source_neutral_Xn_alloc(
            get_idx_range<Species, GridXNeutrals>(neutrals_moments));
    DFieldSpXn particle_source_neutral_on_Xn = get_field(particle_source_neutral_Xn_alloc);
    DFieldSpXn density_neutrals = neutrals_moments[GeometryMX::density_idx];
    get_particle_source_term(
            particle_source_neutral_on_Xn,
            get_const_field(plasma_moments_on_Xn[GeometryMX::density_idx]),
            get_const_field(density_neutrals),
            get_const_field(ionisation_rate_on_Xn),
            get_const_field(recombination_rate_on_Xn));

    // S(v) velocity shape calculation for kinetic species
    DFieldMemSpXVx plasma_source_alloc(get_idx_range(allfdistribu));
    DFieldSpXVx plasma_source = get_field(plasma_source_alloc);
    DFieldSpX plasma_temperature_on_X = plasma_moments_on_X[GeometryMX::temperature_idx];
    get_plasma_source_term(
            plasma_source,
            get_const_field(plasma_temperature_on_X),
            get_const_field(particle_source_neutral_on_Xn));

    // do the actual time stepping
    RK2<DFieldMemSpXVx> timestepper_kinetic(get_idx_range(allfdistribu));
    timestepper_kinetic.update(allfdistribu, dt, [&](DFieldSpXVx df, DConstFieldSpXVx f) {
        get_derivative_allfdistribu(df, get_const_field(plasma_source));
    });

    RK2<DFieldMemSpXn> timestepper_neutrals(get_idx_range(density_neutrals));
    timestepper_neutrals.update(density_neutrals, dt, [&](DFieldSpXn dn, DConstFieldSpXn n) {
        get_derivative_neutrals(dn, get_const_field(particle_source_neutral_on_Xn));
    });

    Kokkos::Profiling::popRegion();
}
