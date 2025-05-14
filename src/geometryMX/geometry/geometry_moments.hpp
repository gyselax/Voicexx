// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "geometry.hpp"

// we want to use GridMom as our first index as we will often pick
// a specific fluid moment and slice at it
using IdxMomSpX = Idx<GridMom, Species, GridX>;
using IdxStepMomSpX = IdxStep<GridMom, Species, GridX>;
using IdxRangeMomSpX = IdxRange<GridMom, Species, GridX>;

template <class ElementType>
using FieldMemMomSpX = FieldMem<ElementType, IdxRangeMomSpX>;
using DFieldMemMomSpX = FieldMemMomSpX<double>;

template <class ElementType>
using FieldMomSpX = Field<ElementType, IdxRangeMomSpX>;

using DFieldMomSpX = FieldMomSpX<double>;

template <class ElementType>
using ConstFieldMomSpX = ConstField<ElementType, IdxRangeMomSpX>;

using DConstFieldMomSpX = ConstFieldMomSpX<double>;

namespace GeometryMX {

/**
 * @brief An enum containing the three values 0, 1 and 2 corresponding to the index
 * of the moments.
 * This is an enum to be able to use it in device code.
 */
enum MomentIdx { density_moment = 0, velocity_moment = 1, temperature_moment = 2 };

static constexpr IdxMom density_idx(density_moment);
static constexpr IdxMom velocity_idx(velocity_moment);
static constexpr IdxMom temperature_idx(temperature_moment);

static constexpr IdxRangeMom first_three_moments(density_idx, IdxStepMom(3));

/**
 * @brief Check whether the index range on the moments given is only containing the density
 * @param[in] moments The index range to check
 */
inline bool is_valid(IdxRangeMom const& moments)
{
    if (moments.size() <= 0 || moments.size() > 3) {
        return false;
    }
    if (moments.front() != density_idx) {
        return false;
    }
    return true;
}

/**
 * @brief Check whether the index range on the moments given is only containing the density
 * @param[in] moments The index range to check
 */
inline bool is_only_density(IdxRangeMom const& moments)
{
    return is_valid(moments) && moments.size() == 1;
}

/**
 * @brief Check whether the index range on the moments given is only containing the density and the velocity
 * @param[in] moments The index range to check
 */
inline bool is_density_and_flux(IdxRangeMom const& moments)
{
    return is_valid(moments) && moments.size() == 2 && moments.back() == velocity_idx;
}
} // namespace GeometryMX
