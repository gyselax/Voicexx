// SPDX-License-Identifier: MIT
#pragma once

#include "geometry.hpp"
#include "geometry_neutrals.hpp"
#include "igridneutralfluidsolver.hpp"

/**
 * @brief A dummy class that solves a fluid model.
 * The fluid model leaves the moments of the fluid species unchanged.
 */
class NullGridNeutralSolver : public IGridNeutralFluidSolver
{
public:
    /**
     * @brief The constructor for the class.
     *
     * @param[in] idx_range_fluidsp The moments index range on which the fluid species is defined.
     */
    explicit NullGridNeutralSolver(IdxRangeSp const& idx_range_fluidsp);

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
    DFieldSpMomXn operator()(
            DFieldSpMomXn fluid_moments,
            DConstFieldSpXVx allfdistribu,
            DConstFieldX efield,
            double dt) const override;
};
