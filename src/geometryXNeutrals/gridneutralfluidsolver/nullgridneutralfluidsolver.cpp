// SPDX-License-Identifier: MIT

#include <stdexcept>

#include "nullgridneutralfluidsolver.hpp"

NullGridNeutralSolver::NullGridNeutralSolver(IdxRangeSp const& idx_range_fluidsp)
{
    // charged fluid species is not allowed for now
    for (IdxSp const isp : idx_range_fluidsp) {
        if (charge(isp) != 0.) {
            throw std::runtime_error("Neutrals charge should be zero");
        }
    }
}

DFieldSpMomXn NullGridNeutralSolver::operator()(
        DFieldSpMomXn const fluid_moments,
        DConstFieldSpXVx const allfdistribu,
        DConstFieldX const efield,
        double const dt) const
{
    return fluid_moments;
}
