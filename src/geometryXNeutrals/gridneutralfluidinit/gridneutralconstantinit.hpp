// SPDX-License-Identifier: MIT

#pragma once

#include "geometry_neutrals.hpp"
#include "igridneutralfluidinit.hpp"
/**
 * @brief A class that initialises on its own grid a fluid species with constant moments.
 */
class GridNeutralConstantInit : public IGridNeutralFluidInit
{
    DFieldMemSpMom m_moments_alloc; // the value of the moments does not depend on x

public:
    /**
     * @brief Creates an instance of the GridNeutralConstantInit class.
     * @param[in] moments The fluid moments the fluid species should be initialised with.
     */
    explicit GridNeutralConstantInit(host_t<DConstFieldSpMom> moments);

    /**
     * @brief Initialises the fluid species with a constant moments.
     * @param[inout] fluid_moments On input: a field referencing an uninitialized fluid species described through its moments.
     *                             On output: a field referencing a the fluid species initialised with constant moments.
     * @return A field referencing the initialised fluid species.
     */
    DFieldSpMomXn operator()(DFieldSpMomXn const fluid_moments) const override;
};
