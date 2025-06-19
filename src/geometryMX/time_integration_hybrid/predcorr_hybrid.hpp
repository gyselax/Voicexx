// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "geometry.hpp"
#include "iboltzmannsolver.hpp"
#include "ifluidsolver.hpp"
#include "iplasmaneutralscoupling.hpp"
#include "iqnsolver.hpp"
#include "itimesolver_hybrid.hpp"

/**
 * @brief A class that solves a Boltzmann-Poisson system of equations coupled to a fluid particles model using a predictor-corrector scheme.
 *
 * A class that solves a Boltzmann-Poisson system with
 * a predictor corrector scheme. This scheme consists in
 * estimating the electric potential after a time interval
 * of a half-timestep. This potential is then used to compute
 * the value of the distribution function at time t+dt, and the
 * value of the fluid moments for the fluid species.
 * The fluid moments value is computed first without the source term S_n,N(x) via the fluid_solver.
 * Then, a coupling between the distribution function and the fluid moments is computed with the kinetic_fluid_coupling function.
 * dt is the timestep of the simulation.
 */
template <typename FluidGrid>
class PredCorrHybrid : public ITimeSolverHybrid<FluidGrid>
{
private:
    using IdxRangeMomSpGrid = IdxRange<GridMom, Species, FluidGrid>;
    using DFieldMomSpGrid = DField<IdxRangeMomSpGrid>;
    using DFieldMemMomSpGrid = DFieldMem<IdxRangeMomSpGrid>;

private:
    IBoltzmannSolver const& m_boltzmann_solver;
    IFluidSolver<FluidGrid> const& m_fluid_solver;
    IQNSolver const& m_poisson_solver;
    IPlasmaNeutralsCoupling<FluidGrid> const& m_kinetic_fluid_coupling;

public:
    /**
     * @brief Creates an instance of the predictor-corrector class.
     * @param[in] boltzmann_solver A solver for a Boltzmann equation.
     * @param[in] fluid_solver A solver for a fluid model.
     * @param[in] poisson_solver A solver for a Quasi-Neutrality equation.
     * @param[in] kinetic_fluid_coupling A solver of the neutral source term in both the Boltzmann and fluid equations.
     */
    PredCorrHybrid(
            IBoltzmannSolver const& boltzmann_solver,
            IFluidSolver<FluidGrid> const& fluid_solver,
            IQNSolver const& poisson_solver,
            IPlasmaNeutralsCoupling<FluidGrid> const& kinetic_fluid_coupling)
        : m_boltzmann_solver(boltzmann_solver)
        , m_fluid_solver(fluid_solver)
        , m_poisson_solver(poisson_solver)
        , m_kinetic_fluid_coupling(kinetic_fluid_coupling)
    {
    }

    /**
     * @brief Solves the Boltzmann-Poisson-fluid system.
     * @param[in, out] allfdistribu On input: the initial value of the distribution function.
     *                              On output: the value of the distribution function after solving
     *                              the Boltzmann-Poisson-fluid system a given number of iterations.
     ** @param[in, out] fluid_moments On input: a field referencing the fluid species.
     *                                On output: the state of the fluid species after solving
     *                                the Boltzmann-Poisson-fluid system a given number of iterations.
     * @param[in] time_start The physical time at the start of the simulation.
     * @param[in] dt The timestep.
     * @param[in] steps The number of iterations to be performed by the predictor-corrector.
     * @return The distribution function after solving the system.
     */
    DFieldSpXVx operator()(
            DFieldSpXVx allfdistribu,
            DFieldMomSpGrid fluid_moments,
            double time_start,
            double dt,
            int steps = 1) const override
    {
        auto allfdistribu_alloc = ddc::create_mirror_view(allfdistribu);
        host_t<DFieldSpXVx> allfdistribu_host = get_field(allfdistribu_alloc);

        IdxRangeX const idx_range_x = get_idx_range<GridX>(allfdistribu);

        // electrostatic potential and electric field (depending only on x)
        host_t<DFieldMemX> electrostatic_potential_host(idx_range_x);
        DFieldMemX electrostatic_potential(idx_range_x);

        DFieldMemX electric_field(idx_range_x);

        host_t<DFieldMemMomSpGrid> fluid_moments_host(get_idx_range(fluid_moments));

        // a 2D chunk of the same size as fdistribu
        host_t<DFieldMemSpXVx> allfdistribu_half_t_host(get_idx_range(allfdistribu));
        DFieldMemSpXVx allfdistribu_half_t(get_idx_range(allfdistribu));

        m_poisson_solver(
                get_field(electrostatic_potential),
                get_field(electric_field),
                get_const_field(allfdistribu));

        int iter = 0;
        for (; iter < steps; ++iter) {
            double const iter_time = time_start + iter * dt;

            // computation of the electrostatic potential at time tn and
            // the associated electric field
            m_poisson_solver(
                    get_field(electrostatic_potential),
                    get_field(electric_field),
                    get_const_field(allfdistribu));
            // copies necessary to PDI
            ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
            ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);
            ddc::parallel_deepcopy(fluid_moments_host, fluid_moments);
            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .with("time_saved", iter_time)
                    .with("fdistribu", allfdistribu_host)
                    .with("fluid_moments", fluid_moments_host)
                    .with("electrostatic_potential", electrostatic_potential_host);

            // copy fdistribu
            ddc::parallel_deepcopy(allfdistribu_half_t, allfdistribu);

            // predictor
            m_boltzmann_solver(
                    get_field(allfdistribu_half_t),
                    get_const_field(electric_field),
                    dt / 2);

            // computation of the electrostatic potential at time tn+1/2
            // and the associated electric field
            m_poisson_solver(
                    get_field(electrostatic_potential),
                    get_field(electric_field),
                    get_const_field(allfdistribu_half_t));
            // correction on a dt
            m_boltzmann_solver(allfdistribu, get_const_field(electric_field), dt);
            m_fluid_solver(
                    fluid_moments,
                    get_const_field(allfdistribu),
                    get_const_field(electric_field),
                    dt);

            m_kinetic_fluid_coupling(allfdistribu, fluid_moments, dt);
        }

        double const final_time = time_start + iter * dt;
        m_poisson_solver(
                get_field(electrostatic_potential),
                get_field(electric_field),
                get_const_field(allfdistribu));

        ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
        ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);
        ddc::parallel_deepcopy(fluid_moments_host, fluid_moments);
        ddc::PdiEvent("last_iteration")
                .with("iter", iter)
                .with("time_saved", final_time)
                .with("fdistribu", allfdistribu_host)
                .with("fluid_moments", fluid_moments_host)
                .with("electrostatic_potential", electrostatic_potential_host);

        return allfdistribu;
    }
};
