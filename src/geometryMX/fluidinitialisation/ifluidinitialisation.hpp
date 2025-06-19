// SPDX-License-Identifier: MIT
#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class that allows for initialising a fluid species.
 * The template is used to precise the grid on which the species should be initialised.
 */
template <typename FluidGrid>
class IFluidInitialisation
{
private:
    using IdxRangeMomSpGrid = IdxRange<GridMom, Species, FluidGrid>;
    using DFieldMomSpGrid = Field<double, IdxRangeMomSpGrid>;

public:
    virtual ~IFluidInitialisation() = default;

    /**
     * @brief Operator for initialising a neutral species.
     * @param[in, out] fluid_moments On input: the uninitialized fluid species.
     *                                 On output: the initialised fluid species.
     * @return A field referencing the initialised fluid species.
     */
    virtual DFieldMomSpGrid operator()(DFieldMomSpGrid fluid_moments) const = 0;
};
