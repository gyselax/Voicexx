# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(sheath_${GEOMETRY_VARIANT} sheath.cpp)
target_link_libraries(sheath_${GEOMETRY_VARIANT}
  PUBLIC
  DDC::core
  DDC::pdi
  gslx::advection
  gslx::boltzmann_${GEOMETRY_VARIANT}
  gslx::initialisation_${GEOMETRY_VARIANT}
  gslx::interpolation
  gslx::io
  gslx::paraconfpp
  gslx::poisson_${GEOMETRY_VARIANT}
  gslx::speciesinfo
  gslx::time_integration_${GEOMETRY_VARIANT}
  gslx::utils
  voicexx::plasmasources_${GEOMETRY_VARIANT}
  voicexx::spline_definitions

  PDI::PDI_C
  paraconf::paraconf
  )
set_target_properties(sheath_${GEOMETRY_VARIANT}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
)

endforeach()
