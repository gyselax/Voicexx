// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "geometry_neutrals.hpp"
#include "gridneutral_interpolator.hpp"
#include "ifluidsolver.hpp"
#include "ireactionrate.hpp"

/**
 * @brief A class that solves a so-called "pressure-diffusive" fluid neutral model.
 *
 * The equation of the model that describes the evolution of the density of
 * neutrals can be written in dimensional units as
 *
 * @f$\partial_t n_n + \partial_x (n_{n,eq} u_i - D_p \partial_x(T_n n_n)) = S_n, @f$
 *
 * where @f$n_n(x,t)@f$ is the time and space dependent neutral density, @f$T_n@f$ is the
 * temperature of neutrals (considered equal to the ion temperature) and @f$u_i(x)@f$ is the ion fluid velocity.
 *
 * In the above equation the following definitions are used:
 *
 * @f$n_{n,eq} = \frac{n_i n_e K_r + n_n n_i K_{cx}}{n_i K_{cx} + n_e K_i}@f$
 *
 * and
 *
 * @f$D_p = \frac{1}{m_n (n_i K_{cx} + n_e K_i)}@f$
 *
 * where @f$n_i@f$ (resp. @f$n_e@f$) is the ion (resp. electron) density and @f$m_n@f$ stands
 * for the mass of neutrals. The @f$K_i@f$, @f$K_r@f$ and @f$K_{cx}@f$ coefficients
 * represent the reaction rates of ionisation, recombination and charge-exchange reactions.
 *
 * The density source term @f$S_n@f$ is chosen to be 0 her, and will be solved
 * by the plasma-neutrals coupling operator.
 *
 *
 * The pressure-diffusive equation is normalised to the relevant normalisation quantities:
 * - densities to a reference density @f$n_0@f$;
 * - temperatures to a reference temperature @f$T_0@f$;
 * - time to the electron plasma frequency @f$\omega_{pe0} = \sqrt{n_0 e^2/(m_e \varepsilon_0)}@f$;
 * - space to the Debye length @f$\lambda_{D0} = \sqrt{\varepsilon_0 T_0 / (n_0 e^2)}@f$;
 * - ion mean velocity to the ion thermal velocity @f$v_{Ti0} = \sqrt{T_0/m_i}@f$;
 * - reaction rates to a reference rate @f$K_{cx,0}@f$;
 * - masses to the electron mass @f$m_e@f$.
 *
 * Moreover, to be able to see an effect of reaction between the plasma and the neutrals,
 * even if the simulation box is way smaller than the reality, the reaction rates are all multiplied
 * by a constant.
 *
 * The pressure-diffusive model is solved using a RK2 time integrator.
 * Spatial derivatives are computed using splines polynomials.
 */
class DiffGridsFluidSolver : public IFluidSolver<GridXNeutrals>
{
public:
    /**
     * @brief An enum class to store the type of boundary condition used
     * It can be:
     * - recycling: the inward neutral flux is equal to the outward ion flux
     * - zero flux: the neutral flux is zero at the boundary
     * - escaping neutrals: no boundary condition for the flux
     */
    enum class NeutralFluxBoundaryCondition { recycling, zero_flux, escaping_neutrals };

    /**
     * @brief A method to convert a string to a valid boundary condition
     * @param[in] bc_flux_input The input string
     * @return The corresponding boundary condition
     */
    static NeutralFluxBoundaryCondition neutral_flux_boundary_condition(
            std::string const& bc_flux_input);

private:
    IReactionRate const& m_charge_exchange;
    IReactionRate const& m_ionisation;
    IReactionRate const& m_recombination;

    double const m_mean_free_path;

    SplineXNeutralsBuilder const& m_spline_builder_on_Xn;
    SplineXn_GridXnEvaluator const& m_spline_evaluator_on_Xn;

    SplineXBuilder const& m_spline_builder_on_X;
    GridNeutralInterpolator const& m_interpolator;

    DConstFieldVx const m_quadrature_coeffs;

    NeutralFluxBoundaryCondition const m_flux_BC;
    double const m_recycling_coefficient;

private:
    IdxSp find_ion(IdxRangeSp const idx_range_kinsp) const;

public:
    /**
     * @brief Creates an instance of the DiffusiveNeutralSolver class.
     * @param[in] charge_exchange An object that represents charge-exchange reaction rate.
     * @param[in] ionisation An object that represents ionisation reaction rate.
     * @param[in] recombination An object that represents recombination reaction rate.
     * @param[in] mean_free_path The mean free path between two charge-exchange reactions.
     * @param[in] spline_builder_on_Xn A one-dimensional spline builder on GridNeutrals
     * @param[in] spline_evaluator_on_Xn A one-dimensional spline evaluator on GridNeutrals
     * @param[in] spline_builder_on_X A one-dimensional spline builder on GridX
     * @param[in] interpolator_from_X_to_Xn An interpolator to communicate between the two grids
     * @param[in] quadrature_coeffs A constant field referencing coefficients for a quadrature.
     * @param[in] flux_BC The boundary condition on the neutrals particle flux.
     * @param[in] recycling_coeff The neutrals recycling coefficient for density.
     */
    DiffGridsFluidSolver(
            IReactionRate const& charge_exchange,
            IReactionRate const& ionisation,
            IReactionRate const& recombination,
            double const mean_free_path,
            SplineXNeutralsBuilder const& spline_builder_on_Xn,
            SplineXn_GridXnEvaluator const& spline_evaluator_on_Xn,
            SplineXBuilder const& spline_builder_on_X,
            GridNeutralInterpolator const& interpolator_from_X_to_Xn,
            DConstFieldVx const& quadrature_coeffs,
            NeutralFluxBoundaryCondition flux_BC,
            double recycling_coeff = 1);

    /**
     * @brief Updates the neutral fluid moments according to the pressure-diffusive neutral model.
     *
     * Within the pressure-diffusive model only the neutral density is evolved thus it is the
     * only fluid moments we consider.
     *
     * @param[inout] neutrals The fluid moments describing the neutrals.
     * @param[in] allfdistribu A constant Field referencing the distribution function.
     * @param[in] efield A constant Field referencing the electric field.
     * @param[in] dt The time step.
     *
     * @return A field referencing the neutral fluid moments passed as argument.
     */
    DFieldSpMomXn operator()(
            DFieldSpMomXn neutrals,
            DConstFieldSpXVx allfdistribu,
            DConstFieldX efield,
            double dt) const override;

    /**
     * @brief Computes the expression of the time derivative of the neutral fluid moments.
     *
     * The expression of the time derivative is given by the equation of the pressure-diffusive
     * neutral model, that is to say
     *
     * @f$\partial_t n_n = - \partial_x (n_{n,eq} u_i - D_p T_n \partial_x n_n)@f$
     *
     * This function is used by the time integrator (RK2 for instance).
     *
     * @param[inout] dn The time derivative of neutral fluid moments.
     * @param[in] n The fluid moments of the neutral species.
     * @param[in] density The plasma density (for ion and electrons).
     * @param[in] velocity The plasma mean velocity (for ion and electrons).
     * @param[in] temperature The plasma temperature (for ion and electrons).
     */
    void get_derivative(
            DFieldSpMomXn dn,
            DConstFieldSpMomXn n,
            DConstFieldSpX density,
            DConstFieldSpX velocity,
            DConstFieldSpX temperature) const;
};
