// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "geometry_neutrals.hpp"
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
private:
    IReactionRate const& m_charge_exchange;
    IReactionRate const& m_ionisation;
    IReactionRate const& m_recombination;

    double const m_mean_free_path;

    SplineXNeutralsBuilder const& m_spline_builder_on_Xn;
    SplineXn_GridXnEvaluator const& m_spline_evaluator_on_Xn;

    SplineXBuilder_1d const& m_spline_builder_on_X;
    SplineX_GridXnEvaluator const& m_interpolator_from_X_to_Xn;

    DConstFieldVx const m_quadrature_coeffs;

    IdxSp find_ion(IdxRangeSp const idx_range_kinsp) const;

public:
    void interpolate_on_neutral_grid(DFieldSpXn field_on_Xn, DConstFieldSpX field_on_X) const;
    /**
     * @brief Creates an instance of the DiffusiveNeutralSolver class.
     * @param[in] charge_exchange An object that represents charge-exchange reaction rate.
     * @param[in] ionisation An object that represents ionisation reaction rate.
     * @param[in] mean_free_path The mean free path between two charge-exchange reactions.
     * @param[in] spline_x_builder_on_Xn A one-dimensional spline builder on GridNeutrals
     * @param[in] spline_x_evaluator_on_Xn A one-dimensional spline evaluator on GridNeutrals
     * @param[in] spline_x_builder_on_X A one-dimensional spline builder on GridX
     * @param[in] interpolator_from_X_to_Xn A one-dimensional spline evaluator with splines on X which can be evaluated on Xn
     * @param[in] quadrature_coeffs A constant field referencing coefficients for a quadrature.
     */
    DiffGridsFluidSolver(
            IReactionRate const& charge_exchange,
            IReactionRate const& ionisation,
            IReactionRate const& recombination,
            double const mean_free_path,
            SplineXNeutralsBuilder const& spline_builder_on_Xn,
            SplineXn_GridXnEvaluator const& spline_evaluator_on_Xn,
            SplineXBuilder_1d const& spline_builder_on_X,
            SplineX_GridXnEvaluator const& intepolator_from_X_to_Xn,
            DConstFieldVx const& quadrature_coeffs);

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
     * @f$\partial_t n_n = - \partial_x (n_{n,eq} u_i - D_p T_n \partial_x n_n)
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
