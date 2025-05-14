// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for solving the transport of a fluid model on a given grid.
 */
template <typename FluidGrid>
class IFluidSolver
{
private:
    using IdxRangeMomSpGrid = IdxRange<GridMom, Species, FluidGrid>;
    using DFieldMomSpGrid = DField<IdxRangeMomSpGrid>;

public:
    /**
     * @brief Operator for solving the fluid model on one timestep.
     * @param[in, out] fluid_moments On input : a field referencing the fluid species.
     *                               On output : a field referencing the fluid species updated
     *                                after solving the fluid model on one timestep.
     * @param[in] allfdistribu A constant field referencing the distribution function.
     * @param[in] efield A constant field referencing the electric field.
     * @param[in] dt The timestep.
     * @return a field referencing the fluid species after solving the fluid model on one timestep.
     */
    virtual DFieldMomSpGrid operator()(
            DFieldMomSpGrid fluid_moments,
            DConstFieldSpXVx allfdistribu,
            DConstFieldX efield,
            double dt) const = 0;

    virtual ~IFluidSolver<FluidGrid>() = default;
};
