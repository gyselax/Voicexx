// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "ifluidinitialisation.hpp"
/**
 * @brief A class that initialises a fluid species with constant moments.
 */
class ConstantFluidInitialisation : public IFluidInitialisation
{
    DFieldMemSpMom m_moments_alloc;

public:
    /**
     * @brief Creates an instance of the ConstantFluidInitialisation class.
     * @param[in] moments The fluid moments the fluid species should be initialised with. 
     */
    ConstantFluidInitialisation(host_t<DConstFieldSpMom> moments);

    ~ConstantFluidInitialisation() override = default;

    /**
     * @brief Initialises the fluid species with a constant moments.
     * @param[inout] fluid_moments On input: a field referencing an uninitialized fluid species described through its moments.
     *                             On output: a field referencing a the fluid species initialised with constant moments.
     * @return A field referencing the initialised fluid species.
     */
    DFieldSpMomX operator()(DFieldSpMomX const fluid_moments) const override;
};
