// SPDX-License-Identifier: MIT
#pragma once

#include "geometry.hpp"
#include "ifluidsolver.hpp"

/**
 * @brief A dummy class that solves a fluid model.
 * The fluid model leaves the moments of the fluid species unchanged.
 */
template <typename FluidGrid>
class NullFluidSolver : public IFluidSolver<FluidGrid>
{
private:
    using IdxRangeSpMomGrid = IdxRange<Species, GridMom, FluidGrid>;
    using DFieldSpMomGrid = DField<IdxRangeSpMomGrid>;

public:
    /**
     * @brief The constructor for the class.
     *
     * @param[in] idx_range_fluidsp The moments index range on which the fluid species is defined.
     */
    explicit NullFluidSolver(IdxRangeSp const& idx_range_fluidsp)
    {
        // charged fluid species is not allowed for now
        for (IdxSp const isp : idx_range_fluidsp) {
            if (charge(isp) != 0.) {
                throw std::runtime_error("Neutrals charge should be zero");
            }
        }
    }

    /**
     * @brief Solves a dummy fluid model on a timestep dt.
     * @param[in, out] fluid_moments On input : a field referencing the moments of the fluid species.
     *                               On output : a field referencing the moments of the fluid species
     *                               updated after solving the dummy fluid model.
     * @param[in] allfdistribu A constant field referencing the distribution function.
     * @param[in] efield A constant field referencing the electric field.
     * @param[in] dt The timestep.
     * @return a field referencing the fluid species after solving the dummy fluid model on one timestep.
     */
    DFieldSpMomGrid operator()(
            DFieldSpMomGrid fluid_moments,
            DConstFieldSpXVx allfdistribu,
            DConstFieldX efield,
            double dt) const override
    {
        return fluid_moments;
    }
};
