// SPDX-License-Identifier: MIT

#include <cassert>

#include "ddc_helper.hpp"
#include "geometry.hpp"
#include "geometry_neutrals.hpp"
#include "gridneutral_interpolator.hpp"

GridNeutralInterpolator::GridNeutralInterpolator(
        SplineXBuilder const& spline_builder_on_X,
        SplineXNeutralsBuilder const& spline_builder_on_Xn,
        SplineX_GridXnEvaluator const& evaluator_from_X_to_Xn,
        SplineXn_GridXEvaluator const& evaluator_from_Xn_to_X)
    : m_spline_builder_on_X(spline_builder_on_X)
    , m_spline_builder_on_Xn(spline_builder_on_Xn)
    , m_evaluator_from_X_to_Xn(evaluator_from_X_to_Xn)
    , m_evaluator_from_Xn_to_X(evaluator_from_Xn_to_X)
{
}

void GridNeutralInterpolator::operator()(DFieldSpXn field_on_Xn, DConstFieldSpX field_on_X) const
{
    assert(get_idx_range<Species>(field_on_Xn) == get_idx_range<Species>(field_on_X));
    DBSFieldMemX spline_coeff_alloc(get_spline_idx_range(m_spline_builder_on_X));
    DBSFieldX spline_coeff(get_field(spline_coeff_alloc));
    ddc::for_each(get_idx_range<Species>(field_on_X), [&](IdxSp const isp) {
        m_spline_builder_on_X(spline_coeff, get_const_field(field_on_X[isp]));
        m_evaluator_from_X_to_Xn(field_on_Xn[isp], get_const_field(spline_coeff));
    });
}

void GridNeutralInterpolator::operator()(DFieldSpX field_on_X, DConstFieldSpXn field_on_Xn) const
{
    assert(get_idx_range<Species>(field_on_Xn) == get_idx_range<Species>(field_on_X));
    DBSFieldMemXn spline_coeff_alloc(get_spline_idx_range(m_spline_builder_on_Xn));
    DBSFieldXn spline_coeff(get_field(spline_coeff_alloc));
    ddc::for_each(get_idx_range<Species>(field_on_Xn), [&](IdxSp const isp) {
        m_spline_builder_on_Xn(spline_coeff, get_const_field(field_on_Xn[isp]));
        m_evaluator_from_Xn_to_X(field_on_X[isp], get_const_field(spline_coeff));
    });
}
