// SPDX-License-Identifier: MIT

#pragma once

#include "geometry_neutrals.hpp"

/**
 * @brief An abstract class that allows for initializing a fluid species on its own grid.
 */
class IGridNeutralFluidInit
{
public:
    /**
     * @brief Operator for initializing a neutral species.
     * @param[in, out] fluid_moments On input: the uninitialized fluid species.
     *                                 On output: the initialized fluid species.
     * @return A field referencing the initialized fluid species.
     */
    virtual DFieldSpMomXn operator()(DFieldSpMomXn fluid_moments) const = 0;
};
