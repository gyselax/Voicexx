// SPDX-License-Identifier: MIT
#pragma once

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  Nx_spline_cells : int
  Nvx_spline_cells : int
  iter : int
  iter_start : int
  time_saved : double
  nbstep_diag: int
  iter_saved : int
  Lx : double
  breakpoints_x_extents: { type: array, subtype: int64, size: 1 }
  breakpoints_x:
    type: array
    subtype: double
    size: [ '$breakpoints_x_extents[0]' ]
  breakpoints_vx_extents: { type: array, subtype: int64, size: 1 }
  breakpoints_vx:
    type: array
    subtype: double
    size: [ '$breakpoints_vx_extents[0]' ]
  MeshX_extents: { type: array, subtype: int64, size: 1 }
  MeshX:
    type: array
    subtype: double
    size: [ '$MeshX_extents[0]' ]
  MeshVx_extents: { type: array, subtype: int64, size: 1 }
  MeshVx:
    type: array
    subtype: double
    size: [ '$MeshVx_extents[0]' ]
  Nkinspecies: int
  fdistribu_charges_extents : { type: array, subtype: int64, size: 1 }
  fdistribu_charges:
    type: array
    subtype: double
    size: [ '$fdistribu_charges_extents[0]' ]
  fdistribu_masses_extents : { type: array, subtype: int64, size: 1 }
  fdistribu_masses:
    type: array
    subtype: double
    size: [ '$fdistribu_masses_extents[0]' ]
  neutrals_masses_extents: { type: array, subtype: int64, size: 1 }
  neutrals_masses:
    type: array
    subtype: double
    size: [ '$neutrals_masses_extents[0]' ]
  fdistribu_eq_extents : { type: array, subtype: int64, size: 2 }
  fdistribu_eq:
    type: array
    subtype: double
    size: [ '$fdistribu_eq_extents[0]', '$fdistribu_eq_extents[1]' ]
  collintra_nustar0 : double
  collinter_nustar0 : double
  normalization_coeff_neutrals : double
  norm_coeff_rate_neutrals : double
  charge_exchange_coefficients:
    type: array
    subtype: double
    size: [ 5 ]
  ionization_slope_coefficients:
    type: array
    subtype: double
    size: [ 6 ]
  ionization_intercept_coefficients:
    type: array
    subtype: double
    size: [ 6 ]
  recombination_slope_coefficients:
    type: array
    subtype: double
    size: [ 2 ]
  recombination_intercept_coefficients:
    type: array
    subtype: double
    size: [ 2 ]

  krook_sink_adaptive_extent : double
  krook_sink_adaptive_stiffness : double
  krook_sink_adaptive_amplitude : double
  krook_sink_adaptive_density : double
  krook_sink_adaptive_temperature : double
  krook_sink_adaptive_mask_extents: { type: array, subtype: int64, size: 1 }
  krook_sink_adaptive_mask:
    type: array
    subtype: double
    size: [ '$krook_sink_adaptive_mask_extents[0]' ]
  krook_sink_adaptive_ftarget_extents: { type: array, subtype: int64, size: 1 }
  krook_sink_adaptive_ftarget:
    type: array
    subtype: double
    size: [ '$krook_sink_adaptive_ftarget_extents[0]' ]
  
  krook_sink_constant_extent : double
  krook_sink_constant_stiffness : double
  krook_sink_constant_amplitude : double
  krook_sink_constant_density : double
  krook_sink_constant_temperature : double
  krook_sink_constant_mask_extents: { type: array, subtype: int64, size: 1 }
  krook_sink_constant_mask:
    type: array
    subtype: double
    size: [ '$krook_sink_constant_mask_extents[0]' ]
  krook_sink_constant_ftarget_extents: { type: array, subtype: int64, size: 1 }
  krook_sink_constant_ftarget:
    type: array
    subtype: double
    size: [ '$krook_sink_constant_ftarget_extents[0]' ]

  kinetic_source_extent : double
  kinetic_source_stiffness : double
  kinetic_source_amplitude : double
  kinetic_source_density : double
  kinetic_source_energy : double
  kinetic_source_temperature : double
  kinetic_source_velocity_shape_extents: { type: array, subtype: int64, size: 1 }
  kinetic_source_velocity_shape:
    type: array
    subtype: double
    size: [ '$kinetic_source_velocity_shape_extents[0]' ]
  kinetic_source_spatial_extent_extents: { type: array, subtype: int64, size: 1 }
  kinetic_source_spatial_extent:
    type: array
    subtype: double
    size: [ '$kinetic_source_spatial_extent_extents[0]' ]

  krook_neutrals_amplitude : double
  krook_neutrals_mask_extents: { type: array, subtype: int64, size: 1 }
  krook_neutrals_mask:
    type: array
    subtype: double
    size: [ '$krook_neutrals_mask_extents[0]' ]

  filename_size: size_t
  filename: {type: array, subtype: char, size: "$filename_size"}

data:
  fdistribu_extents: { type: array, subtype: int64, size: 3 }
  fdistribu:
    type: array
    subtype: double
    size: [ '$fdistribu_extents[0]', '$fdistribu_extents[1]', '$fdistribu_extents[2]' ]
  fluid_moments_extents : { type: array, subtype: int64, size: 3 }
  fluid_moments:
    type: array
    subtype: double
    size: [ '$fluid_moments_extents[0]', '$fluid_moments_extents[1]', '$fluid_moments_extents[2]' ]
  electrostatic_potential_extents: { type: array, subtype: int64, size: 1 }
  electrostatic_potential:
    type: array
    subtype: double
    size: [ '$electrostatic_potential_extents[0]' ]
  charge_exchange_rate_extents: { type: array, subtype: int64, size: 2 }
  charge_exchange_rate:
    type: array
    subtype: double
    size: [ '$charge_exchange_rate_extents[0]' , '$charge_exchange_rate_extents[1]' ]
  ionization_rate_extents: { type: array, subtype: int64, size: 2 }
  ionization_rate:
    type: array
    subtype: double
    size: [ '$ionization_rate_extents[0]' , '$ionization_rate_extents[1]' ]
  recombination_rate_extents: { type: array, subtype: int64, size: 2 }
  recombination_rate:
    type: array
    subtype: double
    size: [ '$recombination_rate_extents[0]' , '$recombination_rate_extents[1]' ]
  diffusion_term_extents: { type: array, subtype: int64, size: 2 }
  diffusion_term:
    type: array
    subtype: double
    size: [ '$diffusion_term_extents[0]' , '$diffusion_term_extents[1]' ]
  convection_term_extents: { type: array, subtype: int64, size: 2 }
  convection_term:
    type: array
    subtype: double
    size: [ '$convection_term_extents[0]' , '$convection_term_extents[1]' ]

plugins:
  set_value:
    on_init:
      - share:
        - iter_saved: 0
    on_data:
      iter:
        - set:
          - iter_saved: '${iter_start} + ${iter}/${nbstep_diag}'
    on_finalize:
      - release: [iter_saved]
  decl_hdf5:
    - file: 'VOICEXX_initstate.h5'
      on_event: [cx_rate_coeff_pol_expose]
      collision_policy: replace_and_warn
      write:
        - charge_exchange_coefficients
    - file: 'VOICEXX_initstate.h5'
      on_event: [i_rate_coeff_pol_expose]
      collision_policy: write_into
      write:
        - ionization_slope_coefficients
        - ionization_intercept_coefficients
    - file: 'VOICEXX_initstate.h5'
      on_event: [r_rate_coeff_pol_expose]
      collision_policy: write_into
      write:
        - recombination_slope_coefficients
        - recombination_intercept_coefficients

    - file: 'VOICEXX_initstate.h5'
      on_event: [initial_state]
      collision_policy: write_into
      write:
        - Nx_spline_cells
        - Nvx_spline_cells
        - MeshX
        - MeshVx
        - Lx
        - nbstep_diag
        - collintra_nustar0
        - collinter_nustar0
        - normalization_coeff_neutrals
        - norm_coeff_rate_neutrals

        - Nkinspecies
        - fdistribu_charges
        - fdistribu_masses
        - neutrals_masses
        - fdistribu_eq

        - krook_sink_adaptive_extent
        - krook_sink_adaptive_stiffness
        - krook_sink_adaptive_amplitude
        - krook_sink_adaptive_density
        - krook_sink_adaptive_temperature
        - krook_sink_adaptive_mask
        - krook_sink_adaptive_ftarget

        - krook_sink_constant_extent
        - krook_sink_constant_stiffness
        - krook_sink_constant_amplitude
        - krook_sink_constant_density
        - krook_sink_constant_temperature
        - krook_sink_constant_mask
        - krook_sink_constant_ftarget

        - kinetic_source_extent
        - kinetic_source_stiffness
        - kinetic_source_amplitude
        - kinetic_source_density
        - kinetic_source_energy
        - kinetic_source_temperature
        - kinetic_source_velocity_shape
        - kinetic_source_spatial_extent

        - krook_neutrals_amplitude
        - krook_neutrals_mask


    - file: 'VOICEXX_${iter_saved:05}.h5'
      on_event: [iteration, last_iteration]
      when: '${iter} % ${nbstep_diag} = 0'
      collision_policy: replace_and_warn
      write: [time_saved, fdistribu, fluid_moments, electrostatic_potential]
    - file: 'VOICEXX_${iter_saved:05}.h5'
      on_event: [reaction_rate_expose]
      when: '${iter} % ${nbstep_diag} = 0'
      collision_policy: write_into
      write: [charge_exchange_rate, ionization_rate, recombination_rate]
    - file: 'VOICEXX_${iter_saved:05}.h5'
      on_event: [diff_conv_expose]
      when: '${iter} % ${nbstep_diag} = 0'
      collision_policy: write_into
      write: [diffusion_term, convection_term]
    - file: 'VOICEXX_${iter_start:05}.h5'
      on_event: restart
      read: [time_saved, fdistribu, fluid_moments]
    - file: '${filename}'
      on_event: [read_x_extents]
      read:
        breakpoints_x_extents: {size_of: breakpoints_x}
    - file: '${filename}'
      on_event: [read_x]
      read:
        breakpoints_x: ~
    - file: '${filename}'
      on_event: [read_vx_extents]
      read:
        breakpoints_vx_extents: {size_of: breakpoints_vx}
    - file: '${filename}'
      on_event: [read_vx]
      read:
        breakpoints_vx: ~
  #trace: ~
)PDI_CFG";
