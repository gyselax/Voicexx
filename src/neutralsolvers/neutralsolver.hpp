// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "geometry.hpp"
#include "ifluidsolver.hpp"
#include "ineutralsolver.hpp"
#include "iplasmaneutralscoupling.hpp"

/**
 * @brief A class that solves the evolution of the neutrals and the impact on the plasma
 * on a given fluid grid.
 * It is composed of two things:
 * a transport solver for the evolution of the neutrals,
 * a plasma-neutrals coupling to compute the source terms for the neutrals and the plasma.
 */
template <typename FluidGrid>
class NeutralSolver : public INeutralSolver<FluidGrid>
{
private:
    using IdxRangeSpMomGrid = IdxRange<Species, GridMom, FluidGrid>;
    using DFieldSpMomGrid = DField<IdxRangeSpMomGrid>;

private:
    IFluidSolver<FluidGrid> const& m_transport_solver;
    IPlasmaNeutralsCoupling<FluidGrid> const& m_plasma_neutrals_coupling;

public:
    /**
     * @brief Creates an instance of the NeutralSolver class.
     * @param[in] transport_solver The solver for the neutrals transport
     * @param[in] plasma_neutrals_coupling The solver for the source terms
     */
    NeutralSolver(
            IFluidSolver<FluidGrid> const& transport_solver,
            IPlasmaNeutralsCoupling<FluidGrid> const& plasma_neutrals_coupling)
        : m_transport_solver(transport_solver)
        , m_plasma_neutrals_coupling(plasma_neutrals_coupling)
    {
    }

    /**
     * @brief Updates the neutral fluid moments according to the transport solver.
     * Updates the neutral fluid moments and the plasma distribution function according
     * to the sources in density, momemtum and energy.
     *
     * @param[inout] neutrals The fluid moments describing the neutrals.
     * @param[inout] allfdistribu A field referencing the distribution function for the plasma.
     * @param[in] dt The time step.
     *
     * @return A field referencing the neutral fluid moments passed as argument.
     */
    DFieldSpMomGrid operator()(DFieldSpMomGrid neutrals, DFieldSpXVx allfdistribu, double dt)
            const override
    {
        m_transport_solver(neutrals, get_const_field(allfdistribu), dt);
        m_plasma_neutrals_coupling(allfdistribu, neutrals, dt);
        return neutrals;
    }
};
