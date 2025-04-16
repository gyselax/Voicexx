#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "densitycoupling.hpp"
#include "geometry.hpp"
#include "geometry_neutrals.hpp"
#include "ireactionrate.hpp"
#include "rk2.hpp"
#include "species_info.hpp"

DensityCoupling::DensityCoupling(
        double const density_coupling_coeff,
        double const momentum_coupling_coeff,
        double const energy_coupling_coeff,
        IReactionRate const& ionization,
        IReactionRate const& recombination,
        SplineXBuilder_1d const& spline_builder_on_X,
        SplineX_GridXnEvaluator const& interpolator_from_X_to_Xn,
        double const mean_free_path,
        DConstFieldVx const& quadrature_coeffs)
    : m_density_coupling_coeff(density_coupling_coeff)
    , m_momentum_coupling_coeff(momentum_coupling_coeff)
    , m_energy_coupling_coeff(energy_coupling_coeff)
    , m_ionization(ionization)
    , m_recombination(recombination)
    , m_mean_free_path(mean_free_path)
    , m_quadrature_coeffs(quadrature_coeffs)
    , m_spline_builder_on_X(spline_builder_on_X)
    , m_interpolator_from_X_to_Xn(interpolator_from_X_to_Xn)
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

void DensityCoupling::get_source_term(
        DFieldSpXn density_source_neutral,
        DConstFieldSpX kinsp_density,
        DConstFieldSpMomXn neutrals,
        DConstFieldSpX ionization,
        DConstFieldSpX recombination) const
{
    // interpolate on the neutral grid
    IdxRangeSpXn idx_range_plasma_onXn(
            get_idx_range<Species>(kinsp_density),
            get_idx_range<GridXNeutrals>(neutrals));
    DFieldMemSpXn kinsp_density_onXn_alloc(idx_range_plasma_onXn);
    DFieldSpXn kinsp_density_onXn = get_field(kinsp_density_onXn_alloc);
    IdxRangeSpXn idx_range_rate_onXn(
            get_idx_range<Species>(ionization),
            get_idx_range<GridXNeutrals>(neutrals));
    DFieldMemSpXn ionization_onXn_alloc(idx_range_rate_onXn);
    DFieldSpXn ionization_onXn = get_field(ionization_onXn_alloc);
    DFieldMemSpXn recombination_onXn_alloc(idx_range_rate_onXn);
    DFieldSpXn recombination_onXn = get_field(recombination_onXn_alloc);
    interpolate_on_neutral_grid(kinsp_density_onXn, kinsp_density);
    interpolate_on_neutral_grid(ionization_onXn, ionization);
    interpolate_on_neutral_grid(recombination_onXn, recombination);

    IdxSp const iion(find_ion(get_idx_range<Species>(kinsp_density)));
    IdxMom idensity(0);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(density_source_neutral),
            KOKKOS_LAMBDA(IdxSpXn const ispxn) {
                IdxXn ixn(ispxn);
                density_source_neutral(ispxn)
                        = -neutrals(idensity, ispxn) * kinsp_density_onXn(ielec(), ixn)
                                  * ionization_onXn(ispxn)
                          + kinsp_density_onXn(iion, ixn) * kinsp_density_onXn(ielec(), ixn)
                                    * recombination_onXn(ispxn);
            });
}
void DensityCoupling::get_plasma_source(
        DFieldSpXVx plasma_source,
        DConstFieldSpX kinsp_temperature,
        DConstFieldSpXn neutral_density_source_on_Xn) const
{
    /*IdxSp const iion(find_ion(get_idx_range<Species>(plasma_source)));*/
    /*double density_coupling_coeff_proxy = m_density_coupling_coeff;*/
    /*double momentum_coupling_coeff_proxy = m_momentum_coupling_coeff;*/
    /*double energy_coupling_coeff_proxy = m_energy_coupling_coeff;*/
    /*double mean_free_path_proxy = m_mean_free_path;*/
    /**/
    // TODO: interpolate and compute the term
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(plasma_source),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                /*            IdxX const ix(ispxvx);*/
                /*            IdxVx const ivx(ispxvx);*/
                /*            CoordVx const coordvx = ddc::coordinate(ivx);*/
                /*            double const neutral_temperature = kinsp_temperature(iion, ix);*/
                /*            double const coordvx_sq = coordvx * coordvx;*/
                /*            double const density_source*/
                /*                    = density_coupling_coeff_proxy*/
                /*                      * (1.5 - coordvx_sq / (2 * neutral_temperature))*/
                /*                      * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));*/
                /*            double const momentum_source*/
                /*                    = momentum_coupling_coeff_proxy * Kokkos::sqrt(2 / neutral_temperature)*/
                /*                      * coordvx * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));*/
                /*            double const energy_source = 2 * energy_coupling_coeff_proxy*/
                /*                                         * (-1 + coordvx_sq / neutral_temperature)*/
                /*                                         * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));*/
                /*            plasma_source(ispxvx) = -(density_source_neutral_on_Xn(ix)*/
                /*                                      / (Kokkos::sqrt(2 * M_PI * neutral_temperature)*/
                /*                                         * normalization_coeff_alpha0_proxy))*/
                /*                                            * density_source*/
                /*                                    + momentum_source + energy_source;*/
                plasma_source(ispxvx) = 0;
            });
}

void DensityCoupling::get_derivative_neutrals(
        DFieldSpMomXn dn,
        DConstFieldSpMomXn neutrals,
        DConstFieldSpXn density_source_neutral,
        double const sqrt_mass_ratio) const
{
    IdxRangeSpXn range_neutrals_spxn(get_idx_range<Species, GridXNeutrals>(neutrals));
    IdxMom const ineutral_density(0);
    double mean_free_path_proxy = m_mean_free_path;
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            range_neutrals_spxn,
            KOKKOS_LAMBDA(IdxSpXn const ifspx) {
                dn(ifspx, ineutral_density)
                        = sqrt_mass_ratio * density_source_neutral(ifspx) / mean_free_path_proxy;
            });
}

void DensityCoupling::get_derivative_allfdistribu(
        DFieldSpXVx df,
        DConstFieldSpXVx allfdistribu,
        DConstFieldSpXVx velocity_shape_source) const
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(allfdistribu),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) { df(ispxvx) = velocity_shape_source(ispxvx); });
}

void DensityCoupling::operator()(
        DFieldSpXVx const allfdistribu,
        DFieldSpMomXn neutrals,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("KineticFluidCouplingSource");

    // kinetic species fluid moments computation
    IdxRangeSpX dom_kspx(get_idx_range(allfdistribu));
    DFieldMemSpX kinsp_density_alloc(dom_kspx);
    DFieldMemSpX kinsp_velocity_alloc(dom_kspx);
    DFieldMemSpX kinsp_temperature_alloc(dom_kspx);

    DFieldSpX kinsp_density = get_field(kinsp_density_alloc);
    DFieldSpX kinsp_velocity = get_field(kinsp_velocity_alloc);
    DFieldSpX kinsp_temperature = get_field(kinsp_temperature_alloc);

    DConstFieldVx quadrature_coeffs = get_const_field(m_quadrature_coeffs);

    IdxRangeVx const idx_range_vx(get_idx_range<GridVx>(allfdistribu));

    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), kinsp_density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_kspx,
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IdxVx const ivx : idx_range_vx) {
                    CoordVx const coordv = ddc::coordinate(ivx);
                    double const val(quadrature_coeffs(ivx) * allfdistribu(ispx, ivx));
                    kinsp_density(ispx) += val;
                    particle_flux += val * coordv;
                    momentum_flux += val * coordv * coordv;
                }
                kinsp_velocity(ispx) = particle_flux / kinsp_density(ispx);
                kinsp_temperature(ispx) = (momentum_flux - particle_flux * kinsp_velocity(ispx))
                                          / kinsp_density(ispx);
            });

    // storing the index ranges
    IdxRangeX grid_x(get_idx_range<GridX>(allfdistribu));
    IdxRangeSpX range_rate_spx(get_idx_range<Species>(neutrals), grid_x);

    // building reaction rates
    DFieldMemSpX ionization_rate_alloc(range_rate_spx);
    DFieldMemSpX recombination_rate_alloc(range_rate_spx);
    DFieldSpX ionization_rate = get_field(ionization_rate_alloc);
    DFieldSpX recombination_rate = get_field(recombination_rate_alloc);
    m_ionization(
            ionization_rate,
            get_const_field(kinsp_density),
            get_const_field(kinsp_temperature));
    m_recombination(
            recombination_rate,
            get_const_field(kinsp_density),
            get_const_field(kinsp_temperature));

    // source term computation, on Xn
    DFieldMemSpXn density_source_neutral_Xn_alloc(get_idx_range<Species, GridXNeutrals>(neutrals));
    DFieldSpXn density_source_neutral_on_Xn = get_field(density_source_neutral_Xn_alloc);
    get_source_term(
            density_source_neutral_on_Xn,
            get_const_field(kinsp_density),
            get_const_field(neutrals),
            get_const_field(ionization_rate),
            get_const_field(recombination_rate));

    // S(v) velocity shape calculation for kinetic species
    DFieldMemSpXVx plasma_source_alloc(get_idx_range(allfdistribu));
    DFieldSpXVx plasma_source = get_field(plasma_source_alloc);
    get_plasma_source(
            plasma_source,
            get_const_field(kinsp_temperature),
            get_const_field(density_source_neutral_on_Xn));

    /*IdxSp const iion(find_ion(get_idx_range<Species>(allfdistribu)));*/
    /*double density_coupling_coeff_proxy = m_density_coupling_coeff;*/
    /*double momentum_coupling_coeff_proxy = m_momentum_coupling_coeff;*/
    /*double energy_coupling_coeff_proxy = m_energy_coupling_coeff;*/
    /*double mean_free_path_proxy=m_mean_free_path;*/
    /**/
    /*ddc::parallel_for_each(*/
    /*        Kokkos::DefaultExecutionSpace(),*/
    /*        get_idx_range(allfdistribu),*/
    /*        KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {*/
    /*            IdxX const ix(ispxvx);*/
    /*            IdxVx const ivx(ispxvx);*/
    /*            CoordVx const coordvx = ddc::coordinate(ivx);*/
    /*            double const neutral_temperature*/
    /*                    = (kinsp_temperature(iion, ix) + kinsp_temperature(ielec(), ix)) / 2.;*/
    /*            double const coordvx_sq = coordvx * coordvx;*/
    /*            double const density_source*/
    /*                    = density_coupling_coeff_proxy*/
    /*                      * (1.5 - coordvx_sq / (2 * neutral_temperature))*/
    /*                      * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));*/
    /*            double const momentum_source*/
    /*                    = momentum_coupling_coeff_proxy * Kokkos::sqrt(2 / neutral_temperature)*/
    /*                      * coordvx * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));*/
    /*            double const energy_source = 2 * energy_coupling_coeff_proxy*/
    /*                                         * (-1 + coordvx_sq / neutral_temperature)*/
    /*                                         * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));*/
    /*            velocity_shape_source(ispxvx) = (density_source_neutral_on_Xn(ix)*/
    /*                                             / (Kokkos::sqrt(2 * M_PI * neutral_temperature)*/
    /*                                                * normalization_coeff_alpha0_proxy))*/
    /*                                                    * density_source*/
    /*                                            + momentum_source + energy_source;*/
    /*        });*/

    // do the actual time stepping
    RK2<DFieldMemSpMomXn> timestepper_neutrals(get_idx_range(neutrals));
    RK2<DFieldMemSpXVx> timestepper_kinetic(get_idx_range(allfdistribu));
    timestepper_kinetic.update(allfdistribu, dt, [&](DFieldSpXVx df, DConstFieldSpXVx f) {
        get_derivative_allfdistribu(df, f, get_const_field(plasma_source));
    });
    IdxSp const iion(find_ion(get_idx_range<Species>(allfdistribu)));
    double const sqrt_mass_ratio(Kokkos::sqrt(mass(ielec()) / mass(iion)));
    timestepper_neutrals.update(neutrals, dt, [&](DFieldSpMomXn dn, DConstFieldSpMomXn n) {
        get_derivative_neutrals(
                dn,
                n,
                get_const_field(density_source_neutral_on_Xn),
                sqrt_mass_ratio);
    });
    Kokkos::Profiling::popRegion();
}

// interpolate to the neutral grid
void DensityCoupling::interpolate_on_neutral_grid(DFieldSpXn field_on_Xn, DConstFieldSpX field_on_X)
        const
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
