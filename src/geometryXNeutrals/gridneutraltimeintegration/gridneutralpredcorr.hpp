// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "geometry_neutrals.hpp"
#include "iboltzmannsolver.hpp"
#include "igridneutralfluidsolver.hpp"
#include "igridneutraltimesolver.hpp"
#include "iplasmaneutralscoupling.hpp"
#include "iqnsolver.hpp"

/**
 * @brief A class that solves a Boltzmann-Poisson system of equations coupled to a fluid particles model using a predictor-corrector scheme.
 * The plasma and the neutrals are not on the same grid, we call the grid for the neutrals GridNeutral
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
class GridNeutralPredCorr : public IGridNeutralTimeSolver
{
private:
    IBoltzmannSolver const& m_boltzmann_solver;

    IGridNeutralFluidSolver const& m_fluid_solver;

    IQNSolver const& m_poisson_solver;

    IPlasmaNeutralsCoupling<GridXNeutrals> const& m_kinetic_fluid_coupling;

public:
    /**
     * @brief Creates an instance of the predictor-corrector class.
     * @param[in] boltzmann_solver A solver for a Boltzmann equation.
     * @param[in] fluid_solver A solver for a fluid model.
     * @param[in] poisson_solver A solver for a Quasi-Neutrality equation.
     * @param[in] kinetic_fluid_coupling A solver of the neutral source term in both the Boltzmann and fluid equations.
     */
    GridNeutralPredCorr(
            IBoltzmannSolver const& boltzmann_solver,
            IGridNeutralFluidSolver const& fluid_solver,
            IQNSolver const& poisson_solver,
            IPlasmaNeutralsCoupling<GridXNeutrals> const& kinetic_fluid_coupling);

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
            DFieldSpMomXn fluid_moments,
            double time_start,
            double dt,
            int steps = 1) const override;
};
