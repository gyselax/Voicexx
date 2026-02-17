// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "geometry_xvx.hpp"
#include "ifluidinitialisation.hpp"

/**
 * @brief A class that initialises a fluid species with constant moments on a given grid.
 */

template <typename FluidGrid>
class ConstantFluidInitialisation : public IFluidInitialisation<FluidGrid>
{
private:
    using IdxRangeSpMomGrid = IdxRange<Species, GridMom, FluidGrid>;
    using DFieldSpMomGrid = Field<double, IdxRangeSpMomGrid>;
    using IdxSpMomGrid = Idx<Species, GridMom, FluidGrid>;

private:
    DFieldMemSpMom m_moments_alloc;

public:
    /**
     * @brief A broadcast function used to initialise the fluid moments.
     * We only need this because of a nvcc bug.
     */

    class BroadcastFn
    {
        DFieldSpMomGrid m_fluid_moments;

        DConstFieldSpMom m_moments;

    public:
        /** @brief Create an instance of BroadcastFn
         *
         * @param[in] fluid_moments The field to initialise.
         * @param[in] moments The values it should be initialised with.
         */
        BroadcastFn(DFieldSpMomGrid const fluid_moments, DConstFieldSpMom const moments)
            : m_fluid_moments(fluid_moments)
            , m_moments(moments)
        {
        }

        /** @brief Initialise the fluid moments at a given index
         *
         * @param[in] ispmx The index we are looking at
         */
        KOKKOS_FUNCTION void operator()(IdxSpMomGrid const ispmx) const
        {
            IdxSpMom ispm(ispmx);
            m_fluid_moments(ispmx) = m_moments(ispm);
        }
    };

    /**
     * @brief Creates an instance of the ConstantFluidInitialisation class.
     * @param[in] moments The fluid moments the fluid species should be initialised with.
     */
    ConstantFluidInitialisation(host_t<DConstFieldSpMom> moments)
        : m_moments_alloc(get_idx_range(moments))
    {
        ddc::parallel_deepcopy(get_field(m_moments_alloc), moments);
    }

    /**
     * @brief Initialises the fluid species with a constant moments.
     * @param[inout] fluid_moments On input: a field referencing an uninitialized fluid species described through its moments.
     *                             On output: a field referencing the fluid species initialised with constant moments.
     * @return A field referencing the initialised fluid species.
     */
    DFieldSpMomGrid operator()(DFieldSpMomGrid const fluid_moments) const override
    {
        DConstFieldSpMom moments(get_const_field(m_moments_alloc));
        BroadcastFn moments_broadcast(fluid_moments, moments);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(fluid_moments),
                moments_broadcast);
        return fluid_moments;
    }
};
