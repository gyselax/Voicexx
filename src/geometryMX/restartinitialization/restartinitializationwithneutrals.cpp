// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "geometry.hpp"
#include "restartinitializationwithneutrals.hpp"

RestartInitializationWithNeutrals::RestartInitializationWithNeutrals(
        int iter_start,
        double& time_start)
    : m_iter_start(iter_start)
    , m_time_start(time_start)
{
}

DFieldSpXVx RestartInitializationWithNeutrals::operator()(
        DFieldSpXVx const allfdistribu,
        DFieldSpMomX const fluid_moments) const
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
