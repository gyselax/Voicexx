// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

struct GridMom : Moments
{
};

using IdxMom = Idx<GridMom>;
using IdxSpMom = Idx<Species, GridMom>;
using IdxSpMomX = Idx<Species, GridMom, GridX>;

using IdxStepMom = IdxStep<GridMom>;
using IdxStepSpMom = IdxStep<Species, GridMom>;
using IdxStepSpMomX = IdxStep<Species, GridMom, GridX>;

using IdxRangeMom = IdxRange<GridMom>;
using IdxRangeSpMom = IdxRange<Species, GridMom>;
using IdxRangeSpMomX = IdxRange<Species, GridMom, GridX>;

template <class ElementType>
using FieldMemSpMom = FieldMem<ElementType, IdxRangeSpMom>;

template <class ElementType>
using FieldMemSpMomX = FieldMem<ElementType, IdxRangeSpMomX>;

using DFieldMemSpMom = FieldMemSpMom<double>;
using DFieldMemSpMomX = FieldMemSpMomX<double>;

template <class ElementType>
using FieldSpMomX = Field<ElementType, IdxRangeSpMomX>;

template <class ElementType>
using FieldSpMom = Field<ElementType, IdxRangeSpMom>;

using DFieldSpMomX = FieldSpMomX<double>;
using DFieldSpMom = FieldSpMom<double>;

template <class ElementType>
using ConstFieldSpMom = ConstField<ElementType, IdxRangeSpMom>;

template <class ElementType>
using ConstFieldSpMomX = ConstField<ElementType, IdxRangeSpMomX>;

template <class ElementType>
using ConstFieldSpMom = ConstField<ElementType, IdxRangeSpMom>;

using DConstFieldSpMom = ConstFieldSpMom<double>;
using DConstFieldSpMomX = ConstFieldSpMomX<double>;
using DConstFieldSpMom = ConstFieldSpMom<double>;

static constexpr IdxMom density_idx(0);
static constexpr IdxMom momentum_idx(1);
static constexpr IdxMom energy_idx(2);
