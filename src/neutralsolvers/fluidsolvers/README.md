# Fluid solvers

The `fluidsolvers` folder contains code that allows solving for models that describe a species considered as fluid. Such a fluid species is described by moments of the distribution function (density, particle flux, energy...) rather than by a complete distribution function. Equations describing the time evolution of such a fluid moments are called *fluid equations*. Such fluid equations describe the conservation of fluid moments: conservation of density, particle flux, energy, etc.

The currently implemented solvers are :
- NullFluidSolver: a dummy solver that does not change the moments of the fluid species
- SameGridFluidSolver: used to have the fluid evolving on the same grid as the plasma
- DiffGridsFluidSolver: used when the fluid is evolving on its own grid
