// SPDX-License-Identifier: MIT
#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class that allows for initialising a fluid species.
 */
class IFluidInitialisation
{
public:
    virtual ~IFluidInitialisation() = default;

    /**
     * @brief Operator for initialising a neutral species.
     * @param[in, out] fluid_moments On input: the uninitialized fluid species.
     *                                 On output: the initialised fluid species.
     * @return A field referencing the initialised fluid species.
     */
    virtual DFieldSpMomX operator()(DFieldSpMomX fluid_moments) const = 0;
};
