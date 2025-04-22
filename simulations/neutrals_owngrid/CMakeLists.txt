# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(neutrals_owngrid_${GEOMETRY_VARIANT} neutrals_owngrid.cpp)
target_link_libraries(neutrals_owngrid_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::advection
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::moments
        gslx::interpolation
        gslx::paraconfpp
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::io
        gslx::utils
        gslx::rhs_${GEOMETRY_VARIANT}
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::rhs_${GEOMETRY_VARIANT}
        voicexx::fluidinitialisation_${GEOMETRY_VARIANT}
        voicexx::fluidrestart_${GEOMETRY_VARIANT}
        voicexx::fluidsolver_${GEOMETRY_VARIANT}
        voicexx::geometry_neutrals_${GEOMETRY_VARIANT}
        voicexx::reactionrates_${GEOMETRY_VARIANT}
        voicexx::time_integration_hybrid_${GEOMETRY_VARIANT}

        paraconf::paraconf
        PDI::pdi
)

set_target_properties(neutrals_owngrid_${GEOMETRY_VARIANT}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
)

endforeach()
