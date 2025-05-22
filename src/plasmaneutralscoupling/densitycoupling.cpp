// SPDX-License-Identifier: MIT

#include <stdexcept>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "ddc_alias_inline_functions.hpp"
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
        double const temperature_normalisation,
        IReactionRate const& charge_exchange,
        IReactionRate const& ionisation,
        IReactionRate const& recombination,
        GridNeutralInterpolator const& interpolator_between_X_and_Xn,
        double const mean_free_path,
        DConstFieldVx const& quadrature_coeffs)
    : m_density_coupling_coeff(density_coupling_coeff)
    , m_momentum_coupling_coeff(momentum_coupling_coeff)
    , m_energy_coupling_coeff(energy_coupling_coeff)
    , m_T_0(temperature_normalisation)
    , m_charge_exchange(charge_exchange)
    , m_ionisation(ionisation)
    , m_recombination(recombination)
    , m_mean_free_path(mean_free_path)
    , m_quadrature_coeffs(quadrature_coeffs)
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
        DFieldSpXn particle_source_on_Xn,
        DFieldSpX particle_source_on_X,
        DConstFieldSpXn kinsp_density,
        DConstFieldSpXn density_neutrals,
        DConstFieldSpXn ionisation,
        DConstFieldSpXn recombination) const
{
    IdxSp const iion(find_ion(get_idx_range<Species>(kinsp_density)));
    double const sqrt_mass_ratio(Kokkos::sqrt(mass(ielec()) / mass(iion)));
    double const mean_free_path_proxy = m_mean_free_path;
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(particle_source_on_Xn),
            KOKKOS_LAMBDA(IdxSpXn const ispxn) {
                IdxXn ixn(ispxn);
                particle_source_on_Xn(ispxn)
                        = -density_neutrals(ispxn) * kinsp_density(ielec(), ixn) * ionisation(ispxn)
                          + kinsp_density(iion, ixn) * kinsp_density(ielec(), ixn)
                                    * recombination(ispxn);
                particle_source_on_Xn(ispxn) *= sqrt_mass_ratio / mean_free_path_proxy;
            });
    m_interpolator_between_X_and_Xn(particle_source_on_X, get_const_field(particle_source_on_Xn));
}

void DensityCoupling::get_plasma_source_term(
        DFieldSpXVx plasma_source_term,
        DConstFieldMomSpX plasma_moments,
        DConstFieldX neutrals_density,
        DConstFieldX neutrals_particle_flux,
        DConstFieldX neutrals_particle_source,
        DConstFieldX charge_exchange_rate,
        DConstFieldX ionisation_rate,
        DConstFieldX recombination_rate) const
{
    /*double density_coupling_coeff_proxy = m_density_coupling_coeff;*/
    /*double momentum_coupling_coeff_proxy = m_momentum_coupling_coeff;*/
    /*double energy_coupling_coeff_proxy = m_energy_coupling_coeff;*/
    /*double mean_free_path_proxy = m_mean_free_path;*/

    IdxSp const iion(find_ion(get_idx_range<Species>(plasma_source_term)));
    double const mass_ion = mass(iion);
    double const mass_elec = mass(ielec());
    DField<IdxRangeXVx> ions_source_term = plasma_source_term[iion];
    DField<IdxRangeXVx> electrons_source_term = plasma_source_term[ielec()];
    DConstFieldSpX plasma_density = plasma_moments[GeometryMX::density_idx];
    DConstFieldSpX plasma_velocity = plasma_moments[GeometryMX::velocity_idx];
    DConstFieldSpX plasma_temperature = plasma_moments[GeometryMX::temperature_idx];
    double const T_0_proxy = m_T_0;
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range<GridX, GridVx>(plasma_source_term),
            KOKKOS_LAMBDA(IdxXVx const ixvx) {
                IdxX const ix(ixvx);
                IdxVx const ivx(ixvx);
                double const coordvx = ddc::coordinate(ivx);
                double const density_source = neutrals_particle_source(ix); // S_{n,N}
                if (abs(density_source) < 1e-14) {
                    ions_source_term(ixvx) = 0;
                    electrons_source_term(ixvx) = 0;
                } else {
                    double const density_neutrals = neutrals_density(ix);
                    double const density_ions = plasma_density(iion, ix);
                    double const density_electrons = plasma_density(ielec(), ix);
                    double const velocity_ions = plasma_velocity(iion, ix);
                    double const temperature_ions = plasma_temperature(iion, ix);
                    double const temperature_electrons = plasma_temperature(ielec(), ix);
                    double const energy_ions = (mass_ion * velocity_ions * velocity_ions
                                                + density_ions * temperature_ions)
                                               / 2;
                    double const K_cx = charge_exchange_rate(ix);
                    double const K_i = ionisation_rate(ix);
                    double const K_r = recombination_rate(ix);
                    // ion species
                    {
                        double const R_E = 0.25;
                        double const momentum_source
                                = mass_ion
                                  * (neutrals_particle_flux(ix)
                                             * (density_ions * K_i + density_ions * K_cx)
                                     - density_ions * density_electrons * K_r
                                     + density_neutrals * density_ions * K_cx);
                        double const velocity_source = momentum_source / density_source; // U_{N,i}
                        double const energy_source
                                = R_E * energy_ions * density_neutrals * density_electrons * K_i
                                  - energy_ions * density_ions * density_electrons * K_r
                                  - mass_ion * velocity_ions * velocity_ions * density_neutrals
                                            * density_ions * K_cx / 2; // S_{E,i}
                        double const temperature_source = 2 * (energy_source / density_source)
                                                          - mass_ion * velocity_source; // T_{N,i}
                        double const normalisation_term
                                = -density_source
                                  / Kokkos::sqrt(2 * M_PI * temperature_source / mass_ion);
                        ions_source_term(ixvx)
                                = normalisation_term
                                  * Kokkos::exp(
                                          -(mass_ion * Kokkos::pow(coordvx - velocity_source, 2))
                                          / (2 * temperature_source));
                    }
                    // electron species
                    {
                        double const temperature_loss_ionisation
                                = (15 + 170 * Kokkos::exp(-temperature_electrons / 2)) / T_0_proxy;
                        double temperature_loss_recombination
                                = 8 * Kokkos::exp(temperature_electrons / 9);
                        if (temperature_loss_recombination < 250)
                            temperature_loss_recombination = 250;
                        temperature_loss_recombination /= T_0_proxy;
                        double const energy_source
                                = -temperature_loss_ionisation * density_electrons
                                          * density_neutrals * K_i
                                  - temperature_loss_recombination * density_electrons
                                            * density_ions * K_r; // S_{E,e}
                        double const temperature_source
                                = 2 * (energy_source / density_source); // T_{N,e}
                        double const normalisation_term
                                = -density_source
                                  / Kokkos::sqrt(2 * M_PI * temperature_source / mass_elec);
                        electrons_source_term(ixvx) = normalisation_term
                                                      * Kokkos::exp(
                                                              -(mass_elec * Kokkos::pow(coordvx, 2))
                                                              / (2 * temperature_source));
                    }
                }
            });
}

void DensityCoupling::get_derivative_neutrals(
        DFieldSpXn dn,
        DConstFieldSpXn particle_source_neutrals) const
{
    IdxRangeSpXn range_neutrals_spxn(get_idx_range(particle_source_neutrals));
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            range_neutrals_spxn,
            KOKKOS_LAMBDA(IdxSpXn const ifspx) { dn(ifspx) = particle_source_neutrals(ifspx); });
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
        DFieldMomSpX const plasma_moments_on_X,
        DFieldMomSpXn const plasma_moments_on_Xn,
        DConstFieldSpXVx const allfdistribu) const
{
    IdxRangeSpX idx_range_kspx(get_idx_range(allfdistribu));
    DFieldSpX density_on_X = plasma_moments_on_X[GeometryMX::density_idx];
    DFieldSpX velocity_on_X = plasma_moments_on_X[GeometryMX::velocity_idx];
    DFieldSpX temperature_on_X = plasma_moments_on_X[GeometryMX::temperature_idx];

    DConstFieldVx quadrature_coeffs = m_quadrature_coeffs;

    IdxRangeVx const idx_range_vx(get_idx_range<GridVx>(allfdistribu));

    // fluid moments computation
    ddc::parallel_fill(plasma_moments_on_X, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_kspx,
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IdxVx const ivx : idx_range_vx) {
                    CoordVx const coordv = ddc::coordinate(ivx);
                    double const val(quadrature_coeffs(ivx) * allfdistribu(ispx, ivx));
                    density_on_X(ispx) += val;
                    particle_flux += val * coordv;
                    momentum_flux += val * coordv * coordv;
                }
                velocity_on_X(ispx) = particle_flux / density_on_X(ispx);
                temperature_on_X(ispx) = (momentum_flux - particle_flux * velocity_on_X(ispx))
                                         / density_on_X(ispx);
            });
    m_interpolator_between_X_and_Xn(
            plasma_moments_on_Xn[GeometryMX::density_idx],
            get_const_field(density_on_X));
    m_interpolator_between_X_and_Xn(
            plasma_moments_on_Xn[GeometryMX::velocity_idx],
            get_const_field(velocity_on_X));
    m_interpolator_between_X_and_Xn(
            plasma_moments_on_Xn[GeometryMX::temperature_idx],
            get_const_field(temperature_on_X));
}

void DensityCoupling::compute_reaction_rates(
        DFieldSpX const charge_exchange_on_X,
        DFieldSpX const ionisation_on_X,
        DFieldSpX const recombination_on_X,
        DFieldSpXn const charge_exchange_on_Xn,
        DFieldSpXn const ionisation_on_Xn,
        DFieldSpXn const recombination_on_Xn,
        DConstFieldMomSpXn const neutrals_moment_on_Xn,
        DConstFieldMomSpX const plasma_moments_on_X) const
{
    IdxRangeSpX idx_range_rates_on_X(
            get_idx_range<Species>(neutrals_moment_on_Xn),
            get_idx_range<GridX>(plasma_moments_on_X));

    // building reaction rates
    // pay attention, these are normalised to Kcx0=10^-14
    DConstFieldSpX plasma_density = plasma_moments_on_X[GeometryMX::density_idx];
    DConstFieldSpX plasma_temperature = plasma_moments_on_X[GeometryMX::temperature_idx];

    m_charge_exchange(charge_exchange_on_X, plasma_density, plasma_temperature);
    m_ionisation(ionisation_on_X, plasma_density, plasma_temperature);
    m_recombination(recombination_on_X, plasma_density, plasma_temperature);

    m_interpolator_between_X_and_Xn(charge_exchange_on_Xn, get_const_field(charge_exchange_on_X));
    m_interpolator_between_X_and_Xn(ionisation_on_Xn, get_const_field(ionisation_on_X));
    m_interpolator_between_X_and_Xn(recombination_on_Xn, get_const_field(recombination_on_X));
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
    IdxRangeMomSpXn idx_range_kinsp_on_Xn(
            GeometryMX::first_three_moments,
            get_idx_range<Species>(allfdistribu),
            get_idx_range<GridXNeutrals>(neutrals_moments));
    DFieldMemMomSpXn plasma_moments_on_Xn_alloc(idx_range_kinsp_on_Xn);
    DFieldMomSpXn plasma_moments_on_Xn(plasma_moments_on_Xn_alloc);
    compute_plasma_moments(
            plasma_moments_on_X,
            plasma_moments_on_Xn,
            get_const_field(allfdistribu));

    // compute the reaction rates
    IdxRangeSpXn idx_range_neutrals_on_Xn(get_idx_range(neutrals_moments));
    DFieldMemSpXn i_rate_Xn_alloc(idx_range_neutrals_on_Xn);
    DFieldMemSpXn r_rate_Xn_alloc(idx_range_neutrals_on_Xn);
    DFieldMemSpXn cx_rate_Xn_alloc(idx_range_neutrals_on_Xn);
    DFieldSpXn ionisation_rate_on_Xn = get_field(i_rate_Xn_alloc);
    DFieldSpXn recombination_rate_on_Xn = get_field(r_rate_Xn_alloc);
    DFieldSpXn charge_exchange_rate_on_Xn = get_field(cx_rate_Xn_alloc);
    IdxRangeSp neutrals_species(get_idx_range<Species>(neutrals_moments));
    IdxRangeSpX idx_range_neutrals_on_X(neutrals_species, get_idx_range<GridX>(allfdistribu));
    DFieldMemSpX i_rate_X_alloc(idx_range_neutrals_on_X);
    DFieldMemSpX r_rate_X_alloc(idx_range_neutrals_on_X);
    DFieldMemSpX cx_rate_X_alloc(idx_range_neutrals_on_X);
    DFieldSpX ionisation_rate_on_X = get_field(i_rate_X_alloc);
    DFieldSpX recombination_rate_on_X = get_field(r_rate_X_alloc);
    DFieldSpX charge_exchange_rate_on_X = get_field(cx_rate_X_alloc);
    compute_reaction_rates(
            charge_exchange_rate_on_X,
            ionisation_rate_on_X,
            recombination_rate_on_X,
            charge_exchange_rate_on_Xn,
            ionisation_rate_on_Xn,
            recombination_rate_on_Xn,
            get_const_field(neutrals_moments),
            get_const_field(plasma_moments_on_X));

    // compute the particle source term
    DFieldMemSpXn particle_source_on_Xn_alloc(idx_range_neutrals_on_Xn);
    DFieldSpXn particle_source_on_Xn = get_field(particle_source_on_Xn_alloc);
    DFieldMemSpX particle_source_on_X_alloc(idx_range_neutrals_on_X);
    DFieldSpX particle_source_on_X = get_field(particle_source_on_X_alloc);
    DFieldSpXn density_neutrals = neutrals_moments[GeometryMX::density_idx];
    get_particle_source_term(
            particle_source_on_Xn,
            particle_source_on_X,
            get_const_field(plasma_moments_on_Xn[GeometryMX::density_idx]),
            get_const_field(density_neutrals),
            get_const_field(ionisation_rate_on_Xn),
            get_const_field(recombination_rate_on_Xn));

    //interpolate the neutrals fluid moments on gridX
    IdxRangeMomSpX
            idx_range_neutrals_on_MomX(GeometryMX::first_two_moments, idx_range_neutrals_on_X);
    DFieldMemMomSpX neutrals_moments_on_X_alloc(idx_range_neutrals_on_MomX);
    DFieldMomSpX neutrals_moments_on_X = get_field(neutrals_moments_on_X_alloc);
    m_interpolator_between_X_and_Xn(
            neutrals_moments_on_X[GeometryMX::density_idx],
            get_const_field(neutrals_moments[GeometryMX::density_idx]));
    m_interpolator_between_X_and_Xn(
            neutrals_moments_on_X[GeometryMX::velocity_idx],
            get_const_field(neutrals_moments[GeometryMX::velocity_idx]));

    // slicing everything at the first neutrals species
    if (neutrals_species.size() != 1) {
        throw std::runtime_error(
                "For the moments the coupling operator only works for one neutrals species");
    }
    IdxSp ineutrals = neutrals_species.front();
    DConstFieldX plasma_particle_source = get_const_field(particle_source_on_X[ineutrals]);
    DConstFieldX neutrals_density
            = get_const_field(neutrals_moments_on_X[GeometryMX::density_idx][ineutrals]);
    DConstFieldX neutrals_particle_flux
            = get_const_field(neutrals_moments_on_X[GeometryMX::velocity_idx][ineutrals]);
    DConstFieldX charge_exchange = get_const_field(charge_exchange_rate_on_X[ineutrals]);
    DConstFieldX ionisation = get_const_field(ionisation_rate_on_X[ineutrals]);
    DConstFieldX recombination = get_const_field(recombination_rate_on_X[ineutrals]);

    // compute the plasma source term
    DFieldMemSpXVx plasma_source_alloc(get_idx_range(allfdistribu));
    DFieldSpXVx plasma_source = get_field(plasma_source_alloc);
    get_plasma_source_term(
            plasma_source,
            get_const_field(plasma_moments_on_X),
            neutrals_density,
            neutrals_particle_flux,
            plasma_particle_source,
            charge_exchange,
            ionisation,
            recombination);

    // do the actual time stepping
    RK2<DFieldMemSpXVx> timestepper_kinetic(get_idx_range(allfdistribu));
    timestepper_kinetic.update(allfdistribu, dt, [&](DFieldSpXVx df, DConstFieldSpXVx f) {
        get_derivative_allfdistribu(df, get_const_field(plasma_source));
    });

    RK2<DFieldMemSpXn> timestepper_neutrals(get_idx_range(density_neutrals));
    timestepper_neutrals.update(density_neutrals, dt, [&](DFieldSpXn dn, DConstFieldSpXn n) {
        get_derivative_neutrals(dn, get_const_field(particle_source_on_Xn));
    });

    Kokkos::Profiling::popRegion();
}
