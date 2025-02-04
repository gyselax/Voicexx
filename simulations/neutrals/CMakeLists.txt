# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(neutrals_${GEOMETRY_VARIANT} neutrals.cpp)
target_link_libraries(neutrals_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::advection
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::fluidinitialization_${GEOMETRY_VARIANT}
        gslx::fluidsolver_${GEOMETRY_VARIANT}
        gslx::moments
        gslx::initialization_neutrals_${GEOMETRY_VARIANT}
        gslx::interpolation
        gslx::paraconfpp
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::reactionrates_${GEOMETRY_VARIANT}
        gslx::rhs_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::time_integration_hybrid_${GEOMETRY_VARIANT}
        gslx::io
        gslx::utils
        gslx::rhs_${GEOMETRY_VARIANT}

        paraconf::paraconf
        PDI::pdi
)

set_target_properties(neutrals_${GEOMETRY_VARIANT}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
)

endforeach()
