// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "geometry.hpp"
#include "species_info.hpp"

/**
 * @brief A class that initialises the distribution function and the fluid moments from a previous simulation on a given grid.
 *
 * A class that triggers a PDI event to read the values of
 * a distribution function saved in a hdf5 file. These
 * values are copied to the field that represents the
 * distribution function.
 * Note that this class does not heritate from. This is due to the fact iinitialisation only handle one single
 * distribution function. We need two, on different ranges.
 */
template <typename FluidGrid>
class RestartInitialisationWithNeutrals
{
private:
    int m_iter_start; /* iteration number to perform the restart from */
    double& m_time_start; /* corresponding simulation time */

public:
    using IdxRangeSpMomGrid = IdxRange<Species, GridMom, FluidGrid>;
    using DFieldSpMomGrid = Field<double, IdxRangeSpMomGrid>;

    /**
     * @brief Create an initialisation object.
     * @param[in] iter_start An integer representing the number of iteration already performed
     *                       to produce the distribution function used to initialise the current simulation.
     * @param[in] time_start The physical time corresponding to iter_start.
     */
    RestartInitialisationWithNeutrals(int iter_start, double& time_start)
        : m_iter_start(iter_start)
        , m_time_start(time_start)
    {
    }

    ~RestartInitialisationWithNeutrals() = default;

    /**
     * @brief Triggers a PDI event to fill the distribution function with values from a hdf5 file.
     * @param[out] allfdistribu The distribution function initialised with the values
     *                          read from an external file.
     * @param[out] fluidmoments The fluid moments initialised with the values
     *                          read from an external file.
     * @return The initialised distribution function.
     */
    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, DFieldSpMomGrid fluid_moments) const
    {
        auto allfdistribu_host = ddc::create_mirror_view_and_copy(get_field(allfdistribu));
        auto fluid_moments_host = ddc::create_mirror_view_and_copy(get_field(fluid_moments));
        ddc::PdiEvent("restart")
                .with("time_saved", m_time_start)
                .with("fdistribu", allfdistribu_host)
                .with("fluid_moments", fluid_moments_host);
        ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
        ddc::parallel_deepcopy(fluid_moments, fluid_moments_host);
        return allfdistribu;
    }
};
