# Plasma-Neutrals coupling

The `plasmaneutralscoupling` folder contains code that allows for the coupling between the fluid equation of neutrals and the kinetic Boltzmann equation of plasma particles.
This link is made by adding two source terms that introduce or remove particles from both systems in identical and opposite quantities.

The currently implemented couplings are :
- NullPlasmaNeutralsCoupling
- NoEnergyExchangeCoupling (only works if the plasma and the neutrals share the same grid)
