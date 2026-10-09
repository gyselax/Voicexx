# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(sheath_${GEOMETRY_VARIANT} sheath.cpp)
target_link_libraries(sheath_${GEOMETRY_VARIANT}
  PUBLIC
  DDC::core
  DDC::pdi
  gyselalibxx::advection
  gyselalibxx::boltzmann_${GEOMETRY_VARIANT}
  gyselalibxx::initialisation_${GEOMETRY_VARIANT}
  gyselalibxx::interpolation
  gyselalibxx::io
  gyselalibxx::paraconfpp
  gyselalibxx::poisson_${GEOMETRY_VARIANT}
  gyselalibxx::speciesinfo
  gyselalibxx::time_integration_${GEOMETRY_VARIANT}
  gyselalibxx::utils
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
