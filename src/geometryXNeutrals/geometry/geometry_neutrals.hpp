// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

int constexpr BSDegreeXNeutrals = 3;
struct BSplinesXNeutrals : ddc::UniformBSplines<X, BSDegreeXNeutrals>
{
};

auto constexpr SplineXNeutralsBoundary
        = X::PERIODIC ? ddc::BoundCond::PERIODIC : ddc::BoundCond::GREVILLE;
using SplineInterpPointsXNeutrals = ddc::GrevilleInterpolationPoints<
        BSplinesXNeutrals,
        SplineXNeutralsBoundary,
        SplineXNeutralsBoundary>;

struct GridXNeutrals : SplineInterpPointsXNeutrals::interpolation_discrete_dimension_type
{
};

using SplineXNeutralsBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesXNeutrals,
        GridXNeutrals,
        SplineXNeutralsBoundary,
        SplineXNeutralsBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineXn_GridXnEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesXNeutrals,
        GridXNeutrals,
#ifdef PERIODIC_RDIMX
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>
#else
        ddc::ConstantExtrapolationRule<X>,
        ddc::ConstantExtrapolationRule<X>
#endif
        >;

using SplineX_GridXnEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridXNeutrals,
#ifdef PERIODIC_RDIMX
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>
#else
        ddc::ConstantExtrapolationRule<X>,
        ddc::ConstantExtrapolationRule<X>
#endif
        >;

using SplineXn_GridXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesXNeutrals,
        GridX,
#ifdef PERIODIC_RDIMX
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>
#else
        ddc::NullExtrapolationRule,
        ddc::NullExtrapolationRule
#endif
        >;

using IdxXn = Idx<GridXNeutrals>;
using IdxSpXn = Idx<Species, GridXNeutrals>;
using IdxMomXn = Idx<GridMom, GridXNeutrals>;
using IdxSpMomXn = Idx<Species, GridMom, GridXNeutrals>;

using IdxStepXn = IdxStep<GridXNeutrals>;
using IdxStepSpXn = IdxStep<Species, GridXNeutrals>;
using IdxStepMomXn = IdxStep<GridMom, GridXNeutrals>;
using IdxStepSpMomXn = IdxStep<Species, GridMom, GridXNeutrals>;

using IdxRangeBSXn = IdxRange<BSplinesXNeutrals>;

using IdxRangeXn = IdxRange<GridXNeutrals>;
using IdxRangeSpXn = IdxRange<Species, GridXNeutrals>;
using IdxRangeMomXn = IdxRange<GridMom, GridXNeutrals>;
using IdxRangeSpMomXn = IdxRange<Species, GridMom, GridXNeutrals>;

template <class ElementType>
using FieldMemXn = FieldMem<ElementType, IdxRangeXn>;

template <class ElementType>
using BSFieldMemXn = FieldMem<ElementType, IdxRangeBSXn>;

template <class ElementType>
using FieldMemSpXn = FieldMem<ElementType, IdxRangeSpXn>;

template <class ElementType>
using FieldMemSpMomXn = FieldMem<ElementType, IdxRangeSpMomXn>;

using DFieldMemXn = FieldMemXn<double>;
using DBSFieldMemXn = BSFieldMemXn<double>;
using DFieldMemSpXn = FieldMemSpXn<double>;
using DFieldMemSpMomXn = FieldMemSpMomXn<double>;

template <class ElementType>
using BSFieldXn = Field<ElementType, IdxRangeBSXn>;

template <class ElementType>
using FieldXn = Field<ElementType, IdxRangeXn>;

template <class ElementType>
using FieldSpMomXn = Field<ElementType, IdxRangeSpMomXn>;

template <class ElementType>
using FieldSpXn = Field<ElementType, IdxRangeSpXn>;

using DFieldXn = FieldXn<double>;
using DBSFieldXn = BSFieldXn<double>;
using DFieldSpXn = FieldSpXn<double>;
using DFieldSpMomXn = FieldSpMomXn<double>;

template <class ElementType>
using ConstFieldXn = Field<ElementType const, IdxRangeXn>;

template <class ElementType>
using BSConstFieldXn = ConstField<ElementType, IdxRangeBSXn>;

template <class ElementType>
using ConstFieldSpMomXn = ConstField<ElementType, IdxRangeSpMomXn>;

template <class ElementType>
using ConstFieldSpXn = ConstField<ElementType, IdxRangeSpXn>;

using DConstFieldXn = ConstFieldXn<double>;
using DBSConstFieldXn = BSConstFieldXn<double>;
using DConstFieldSpMomXn = ConstFieldSpMomXn<double>;
using DConstFieldSpXn = ConstFieldSpXn<double>;
