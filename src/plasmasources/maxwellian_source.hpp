// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>

#include "geometry.hpp"
#include "irighthandside.hpp"

/**
 * @brief A class that describes a maxwellian source of particles,
 * which can be different for ions and electrons.
 *
 * The MaxwellianSource class solves the following evolution equation:
 * df/dt = S_M
 * Where S_M does not depend on time, and has a maxwellian distribution in velocity.
 * We have f(t+dt) = f(t) + S_kin*dt as a solution of this evolution equation.
 *
 * spatial_extent defines the location where the source is active.
 * spatial_extent is normalised, so that its integral along the
 * spatial direction equals one. It has a hyperbolic tangent shape.
 * It is equal to one in a central zone of the plasma of width defined by the extent parameter.
 */
class MaxwellianSource : public IRightHandSide
{
private:
    double m_amplitude;
    double m_density;
    double m_temperature_elec;
    double m_temperature_ions;
    DFieldMemX m_spatial_extent;
    DFieldMemSpVx m_velocity_shape;

public:
    /**
     * @brief Creates an instance of the MaxwellianSource class.
     * @param[in] idx_range_x The mesh in the x direction.
     * @param[in] idx_range_vx The mesh in the vx direction.
     * @param[in] idx_range_species The range containing all the kinetic species
     * @param[in] extent A parameter that sets the spatial extent of the source.
     * @param[in] stiffness A parameter that sets stiffness of the source extent.
     * @param[in] amplitude A parameter that sets the amplitude of the source.
     * @param[in] density A parameter that sets the density of the source.
     * @param[in] temperature_elec A parameter that sets the electron temperature of the source.
     * @param[in] temperature_ion A parameter that sets the ion temperature of the source.
     */
    MaxwellianSource(
            IdxRangeX const& idx_range_x,
            IdxRangeVx const& idx_range_vx,
            IdxRangeSp const& idx_range_species,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double temperature_elec,
            double temperature_ion);

    /**
     * @brief Update the distribution function following the MaxwellianSource operator.
     *
     * Update the distribution function for both electrons and ions to show how
     * it is modified following the effect of the MaxwellianSource operator.
     *
     * @param[inout] allfdistribu The distribution function.
     * @param[in] dt The time step over which the collisions occur.
     *
     * @return A field referencing the distribution function passed as argument.
     */
    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double dt) const override;
};

namespace maxwellian_source {
MaxwellianSource init_from_input(IdxRangeSpXVx grid_idx_range, PC_tree_t const& yaml_input_file);
}; // namespace maxwellian_source
