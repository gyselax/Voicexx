// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "geometry_neutrals.hpp"

/**
 * @brief A class that is used to interpolate a field from the plasma grid to
 * the neutrals grid or the other way around.
 */

class GridNeutralInterpolator
{
private:
    SplineXBuilder const& m_spline_builder_on_X;
    SplineXNeutralsBuilder const& m_spline_builder_on_Xn;
    SplineX_GridXnEvaluator const& m_evaluator_from_X_to_Xn;
    SplineXn_GridXEvaluator const& m_evaluator_from_Xn_to_X;

public:
    /**
     * @brief Creates an instance of the GridNeutralInterpolator class.
     * @param[in] spline_builder_on_X A spline builder on GridX
     * @param[in] spline_builder_on_Xn A spline builder on GridXNeutrals
     * @param[in] evaluator_from_X_to_Xn A spline evaluator with splines on GridX and evaluation points on GridXNeutrals
     * @param[in] evaluator_from_Xn_to_X A spline evaluator with splines on GridXNeutrals and evaluation points on GridX
     */
    GridNeutralInterpolator(
            SplineXBuilder const& spline_builder_on_X,
            SplineXNeutralsBuilder const& spline_builder_on_Xn,
            SplineX_GridXnEvaluator const& evaluator_from_X_to_Xn,
            SplineXn_GridXEvaluator const& evaluator_from_Xn_to_X);

    /**
     * @brief Do the interpolation from GridX to GridXNeutrals
     * @param[out] field_on_Xn The result of the interpolation on GridXNeutrals
     * @param[in] field_on_X The field on GridX
     */
    void operator()(DFieldSpXn field_on_Xn, DConstFieldSpX field_on_X) const;

    /**
     * @brief Do the interpolation from GridXNeutrals to GridX
     * @param[out] field_on_X The result of the interpolation on GridX
     * @param[in] field_on_Xn The field on GridXNeutrals
     */
    void operator()(DFieldSpX field_on_X, DConstFieldSpXn field_on_Xn) const;
};
