#include "geometry.hpp"
#include "geometry_neutrals.hpp"
#include "nullgridneutralcoupling.hpp"

NullGridNeutralCoupling::NullGridNeutralCoupling() {}

void NullGridNeutralCoupling::operator()(
        DFieldSpXVx const allfdistribu,
        DFieldSpMomXn const neutrals,
        double const dt) const
{
}
