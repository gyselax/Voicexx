// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "geometry_neutrals.hpp"
#include "gridneutral_interpolator.hpp"
#include "iplasmaneutralscoupling.hpp"
#include "ireactionrate.hpp"

/**
 * @brief A class that describes a source of particles due to neutrals.
 *
 * The KineticFluidCouplingSource class solves the following evolution equations:
 * @f$df/dt = S_n,N(x) * S_v(x,v)@f$
 * where @f$S_n,N(x)@f$ is what we call the density_source_neutral
 * @f$dn_N/dt = - S_n,N(x)@f$
 * Where @f$S_n,N(x) = n_N(x) n_e(x) K_i(x) - n_i(x) n_e(x) K_r(x)@f$
 * @f$S_v(x,v)@f$ is the sum of order 0 to 2 Hermite polynomials times a Maxwellian velocity distribution function.
 *
 *
 * The velocity_shape_source @f$S_v(x,v)@f$ defines the velocity profile of the source in the parallel velocity direction.
 * It is the sum of a source that injects only density, a source that injects only momentum and a source that injects only energy.
 * If the density and energy parameters are equal to one (usual case), the resulting velocity_shape is maxwellian.
 *
 * The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/main/doc/geometryXVx/kinetic_source.pdf).
 */
class DensityCoupling : public IPlasmaNeutralsCoupling<GridXNeutrals>
{
private:
    double const m_density_coupling_coeff;
    double const m_momentum_coupling_coeff;
    double const m_energy_coupling_coeff;
    IReactionRate const& m_ionization;
    IReactionRate const& m_recombination;
    double const m_mean_free_path;
    DConstFieldVx const m_quadrature_coeffs;

    SplineXBuilder const& m_spline_builder_on_X;
    GridNeutralInterpolator const& m_interpolator_between_X_and_Xn;

public:
    /**
     * @brief Creates an instance of the KineticFluidCouplingSource class.
     *
     * @param[in] density_coupling_coeff The coefficient of the density source.
     * @param[in] momentum_coupling_coeff The coefficient of the momentum source.
     * @param[in] energy_coupling_coeff The coefficient of the energy source.
     * @param[in] ionization The rate of the ionization reaction.
     * @param[in] recombination The rate of the recombination reaction.
     * @param[in] normalization_coeff The normalization coefficient of neutrals.
     * @param[in] quadrature_coeffs A constant field referencing coefficients for a quadrature.
     * @param[in] mask_extent The extent of the mask for the neutrals fluid.
     * @param[in] mask_stiffnes The stiffnes of the mask for the neutrals fluid.
     * @param[in] gridx The grid on which to construct the wall.
     */
    DensityCoupling(
            double density_coupling_coeff,
            double momentum_coupling_coeff,
            double energy_coupling_coeff,
            IReactionRate const& ionization,
            IReactionRate const& recombination,
            SplineXBuilder const& spline_builder_on_X,
            GridNeutralInterpolator const& interpolator_between_X_and_Xn,
            double const mean_free_path,
            DConstFieldVx const& quadrature_coeffs);

    /**
     * @brief Update the distribution function and neutral density with respect to the density source
     * of each.
     *
     * @param[in, out] allfdistribu The distribution function.
     * @param[in, out] neutrals The neutral density.
     * @param[in] dt The time step over which the collisions occur.
     *
     */
    void operator()(DFieldSpXVx const allfdistribu, DFieldMomSpXn neutrals, double const dt)
            const override;

    /**
     * @brief Computes the source term density_source_neutral(x), with is the result
     * of the sink due to ionization and the source due to recombination
     *
     * @param[in, out] density_source_neutral The source term.
     * @param[in] kinsp_density The computed plasma densities.
     * @param[in] neutrals The neutral density.
     * @param[in] ionization The ionization rate.
     * @param[in] recombination The recombination rate.
     *
    */
    void get_source_term(
            DFieldSpXn density_source_neutral,
            DConstFieldSpX kinsp_density,
            DConstFieldMomSpXn neutrals,
            DConstFieldSpX ionization,
            DConstFieldSpX recombination) const;

    /**
     * @brief the derivative of the neutral density due to the source term
     *
     * @param[in, out] dn The infinitesimal variation of the neutral density.
     * @param[in] neutrals The neutral density.
     * @param[in] density_source_neutral The density source term.
     *
    */
    void get_derivative_neutrals(
            DFieldMomSpXn dn,
            DConstFieldMomSpXn neutrals,
            DConstFieldSpXn density_source_neutral,
            double const sqrt_mass_ratio) const;

    /**
     * @brief Computes df for the equation df/dt = density_source_neutral(x) * velocity_shape_source(x,v).
     *
     * @param[in, out] df The infinitesimal variation of the distribution function.
     * @param[in] allfdistribu The distribution function.
     * @param[in] velocity_shape_source The velocity shape of the source.
     *
    */
    void get_derivative_allfdistribu(
            DFieldSpXVx df,
            DConstFieldSpXVx allfdistribu,
            DConstFieldSpXVx velocity_shape_source) const;

    void get_plasma_source(
            DFieldSpXVx plasma_source,
            DConstFieldSpX kinsp_temperature,
            DConstFieldSpXn neutral_density_source_on_Xn) const;
    void interpolate_on_neutral_grid(DFieldSpXn field_on_Xn, DConstFieldSpX field_on_X) const;

private:
    /**
    * @brief Returns the index of the ion species in the index range.
    *
    * @param[in] idx_range_kinsp The index range of the kinetic species.
    *
    * @return The index of the ion species in the index range.
    */
    IdxSp find_ion(IdxRangeSp const idx_range_kinsp) const;
};
