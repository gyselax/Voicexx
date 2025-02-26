// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "geometry_neutrals.hpp"
#include "gridneutralconstantinit.hpp"

GridNeutralConstantInit::GridNeutralConstantInit(host_t<DConstFieldSpMom> moments)
    : m_moments_alloc(get_idx_range(moments))
{
    ddc::parallel_deepcopy(m_moments_alloc, moments);
}

DFieldSpMomXn GridNeutralConstantInit::operator()(DFieldSpMomXn const fluid_moments) const
{
    DConstFieldSpMom moments(get_field(m_moments_alloc));
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(fluid_moments),
            KOKKOS_LAMBDA(IdxSpMomXn const ispmx) {
                IdxSpMom ispm(ispmx);
                fluid_moments(ispmx) = moments(ispm);
            });
    return fluid_moments;
}
