// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>

#include "geometry.hpp"
#include "iplasmaneutralscoupling.hpp"

/**
 * @brief A class that describes a trivially null source of particles due to neutrals.
 *
 */
template <typename FluidGrid>
class NullPlasmaNeutralsCoupling : public IPlasmaNeutralsCoupling<FluidGrid>
{
private:
    using IdxRangeSpMomGrid = IdxRange<Species, GridMom, FluidGrid>;
    using DFieldSpMomGrid = DField<IdxRangeSpMomGrid>;

public:
    /**
     * @brief Creates an instance of the NullGridNeutralCoupling class.
     *
     */
    NullPlasmaNeutralsCoupling() {}

    /**
     * @brief Update the distribution function and neutral density following the NullGridNeutralCoupling operator.
     * Concretely, in this case, do nothing.
     *
     * @param[in, out] allfdistribu The distribution function.
     * @param[in, out] neutrals The neutral density.
     * @param[in] dt The time step over which the collisions occur.
     *
     */
    void operator()(DFieldSpXVx const allfdistribu, DFieldSpMomGrid const neutrals, double const dt)
            const override
    {
    }
};
