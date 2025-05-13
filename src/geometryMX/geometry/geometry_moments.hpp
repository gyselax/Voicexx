// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "geometry.hpp"

namespace GeometryMX {

static constexpr IdxMom density_idx(0);
static constexpr IdxMom momentum_idx(1);
static constexpr IdxMom energy_idx(2);

static constexpr IdxRangeMom first_three_moments(density_idx, IdxStepMom(3));

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

inline bool is_only_density(IdxRangeMom const& moments)
{
    return is_valid(moments) && moments.size() == 1;
}

inline bool is_density_and_flux(IdxRangeMom const& moments)
{
    return is_valid(moments) && moments.size() == 2 && moments.back() == momentum_idx;
}
} // namespace GeometryMX
