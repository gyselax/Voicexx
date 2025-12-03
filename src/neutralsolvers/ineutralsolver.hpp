// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for solving the transport of neutrals and the interaction between the
 * plasma and the neutrals on a given grid.
 */
template <typename FluidGrid>
class INeutralSolver
{
private:
    using IdxRangeSpMomGrid = IdxRange<Species, GridMom, FluidGrid>;
    using DFieldSpMomGrid = DField<IdxRangeSpMomGrid>;

public:
    virtual ~INeutralSolver<FluidGrid>() = default;

    /**
     * @brief Operator for solving the fluid model on one timestep.
     * @param[in, out] fluid_moments On input : a field referencing the fluid species.
     *                               On output : a field referencing the fluid species updated
     *                                after solving the fluid model on one timestep.
     * @param[in, out] allfdistribu On input: a field referencing the distribution function for the plasma.
     *                              On output: the plasma distribution function updated with the source
     *                              term coming from the neutrals.
     * @param[in] dt The timestep.
     * @return a field referencing the fluid species after solving the fluid model on one timestep.
     */
    virtual DFieldSpMomGrid operator()(
            DFieldSpMomGrid fluid_moments,
            DFieldSpXVx allfdistribu,
            double dt) const = 0;
};
