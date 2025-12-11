// SPDX-License-Identifier: MIT

#pragma once

#ifdef INPUT_MESH
constexpr char const* const mesh_params_yaml = R"PARAMS_CFG(SplineMesh:
  grid_file: "grid.h5"

)PARAMS_CFG";
#else
constexpr char const* const mesh_params_yaml = R"PARAMS_CFG(SplineMesh:
  x_min: 0.0
  x_max: 50
  x_ncells: 512
  vx_min: -6.0
  vx_max: +6.0
  vx_ncells: 256

)PARAMS_CFG";
#endif

constexpr char const* const params_yaml = R"PARAMS_CFG(SpeciesInfo:
- charge: -1.
  mass: 1.
  density_eq: 1.
  temperature_eq: 1.
  mean_velocity_eq: 0.
  perturb_amplitude: 0.
  perturb_mode: 1

- charge: 1.
  mass: 400.
  density_eq: 1.
  temperature_eq: 1.
  mean_velocity_eq: 0.
  perturb_amplitude: 0.
  perturb_mode: 1

NeutralSpeciesInfo:
- mass: 400.
  density_eq: 0.

Krook:
- name: 'constant' # 'constant' or adaptive': constant values or not for nu coeff.
  type: 'sink'
  solver: 'rk2' # possible values : 'rk2'
  extent: 0.10
  stiffness: 1
  amplitude: 0.1
  density: 1e-9
  temperature: 0.5

PlasmaSource:
  extent: 0.45
  stiffness: 4
  amplitude: 0.1
  density: 1.
  temperature_elec: 1.
  temperature_ions: 1.

NeutralsReactionRates:
  temperature_normalisation_in_eV: 10 # in eV
  density_normalisation: 1e20 # n_0 to compute K_{cx,0}
  mean_free_path: 100

TransportSolver:
  on: false

PlasmaNeutralsCoupling:
  on: true
  density_coupling_coeff: 1.0
  momentum_coupling_coeff: 0.0
  energy_coupling_coeff: 0.0

NeutralKrook:
  extent: 0.10
  stiffness: 1.
  amplitude: 0.1

CollisionsInfo:
  enable_inter: true
  nustar0: 0.1

Algorithm:
  deltat: 0.1
  nbiter: 50

Output:
  time_diag: 0.1
)PARAMS_CFG";
