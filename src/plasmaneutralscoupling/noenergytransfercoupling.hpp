// SPDX-License-Identifier: MIT
#pragma once

#include "geometry.hpp"
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
class NoEnergyExchangeCoupling : public IPlasmaNeutralsCoupling<GridX>
{
private:
    double m_density_coupling_coeff;
    double m_momentum_coupling_coeff;
    double m_energy_coupling_coeff;
    IReactionRate const& m_ionisation;
    IReactionRate const& m_recombination;
    double m_normalisation_coeff;
    DConstFieldVx const m_quadrature_coeffs;
    DFieldMemX m_mask;

public:
    /**
     * @brief Creates an instance of the KineticFluidCouplingSource class.
     *
     * @param[in] density_coupling_coeff The coefficient of the density source.
     * @param[in] momentum_coupling_coeff The coefficient of the momentum source.
     * @param[in] energy_coupling_coeff The coefficient of the energy source.
     * @param[in] ionisation The rate of the ionisation reaction.
     * @param[in] recombination The rate of the recombination reaction.
     * @param[in] normalisation_coeff The normalisation coefficient of neutrals.
     * @param[in] quadrature_coeffs A constant field referencing coefficients for a quadrature.
     * @param[in] mask_extent The extent of the mask for the neutrals fluid.
     * @param[in] mask_stiffness The stiffnes of the mask for the neutrals fluid.
     * @param[in] gridx The grid on which to construct the wall.
     */
    NoEnergyExchangeCoupling(
            double density_coupling_coeff,
            double momentum_coupling_coeff,
            double energy_coupling_coeff,
            IReactionRate const& ionisation,
            IReactionRate const& recombination,
            double normalisation_coeff,
            DConstFieldVx const& quadrature_coeffs,
            double mask_extent,
            double mask_stiffness,
            IdxRangeX const& gridx);

    /**
     * @brief Update the distribution function and neutral density following the KineticFluidCouplingSource operator.
     *
     * @param[in, out] allfdistribu The distribution function.
     * @param[in, out] neutrals The neutral density.
     * @param[in] dt The time step over which the collisions occur.
     *
     */
    void operator()(DFieldSpXVx const allfdistribu, DFieldSpMomX neutrals, double const dt)
            const override;

public:
    /**
    * @brief Returns the index of the ion species in the index range.
    *
    * @param[in] idx_range_kinsp The index range of the kinetic species.
    *
    * @return The index of the ion species in the index range.
    */
    IdxSp find_ion(IdxRangeSp const idx_range_kinsp) const;

    /**
     * @brief Computes the source term density_source_neutral(x) of the KineticFluidCouplingSource operator.
     *
     * @param[in, out] density_source_neutral The source term.
     * @param[in] kinsp_density The computed plasma densities of the distribution function.
     * @param[in] neutrals The neutral density.
     * @param[in] ionisation The ionisation rate.
     * @param[in] recombination The recombination rate.
     *
    */
    void get_source_term(
            DFieldX density_source_neutral,
            DConstFieldSpX kinsp_density,
            DConstFieldSpMomX neutrals,
            DConstFieldSpX ionisation,
            DConstFieldSpX recombination) const;

    /**
     * @brief Computes dn (of neutral density) for the equation dn/dt = - density_source_neutral.
     *
     * @param[in, out] dn The infinitesimal of the neutral density.
     * @param[in] neutrals The neutral density.
     * @param[in] density_source_neutral The density source term.
     *
    */
    void get_derivative_neutrals(
            DFieldSpMomX dn,
            DConstFieldSpMomX neutrals,
            DConstFieldX density_source_neutral) const;

    /**
     * @brief Computes df for the equation df/dt = density_source_neutral(x) * velocity_shape_source(x,v).
     *
     * @param[in, out] df The infinitesimal of the distribution function.
     * @param[in] allfdistribu The distribution function.
     * @param[in] velocity_shape_source The velocity shape of the source.
     *
    */
    void get_derivative_allfdistribu(
            DFieldSpXVx df,
            DConstFieldSpXVx allfdistribu,
            DConstFieldSpXVx velocity_shape_source) const;
};
