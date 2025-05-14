// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <stdexcept>

#include "geometry.hpp"
#include "geometry_moments.hpp"
#include "iplasmaneutralscoupling.hpp"

/**
 * @brief A class that describes a trivially null source of particles due to neutrals.
 *
 */
template <typename FluidGrid>
class NullPlasmaNeutralsCoupling : public IPlasmaNeutralsCoupling<FluidGrid>
{
private:
    using IdxRangeMomSpGrid = IdxRange<GridMom, Species, FluidGrid>;
    using DFieldMomSpGrid = DField<IdxRangeMomSpGrid>;

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
    void operator()(DFieldSpXVx const allfdistribu, DFieldMomSpGrid const neutrals, double const dt)
            const override
    {
        if (!GeometryMX::is_valid(ddc::get_domain<GridMom>(neutrals))) {
            throw std::runtime_error("Fluid moments not valid");
        }
    }
};
