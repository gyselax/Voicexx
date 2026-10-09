# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(neutrals_owngrid_${GEOMETRY_VARIANT} neutrals_owngrid.cpp)
target_link_libraries(neutrals_owngrid_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::core
        DDC::pdi
        gyselalibxx::advection
        gyselalibxx::boltzmann_${GEOMETRY_VARIANT}
        gyselalibxx::moments
        gyselalibxx::interpolation
        gyselalibxx::paraconfpp
        gyselalibxx::poisson_${GEOMETRY_VARIANT}
        gyselalibxx::speciesinfo
        gyselalibxx::io
        gyselalibxx::utils
        gyselalibxx::rhs_${GEOMETRY_VARIANT}
        gyselalibxx::geometry_${GEOMETRY_VARIANT}
        gyselalibxx::rhs_${GEOMETRY_VARIANT}
        voicexx::fluidinitialisation_${GEOMETRY_VARIANT}
        voicexx::fluidrestart_${GEOMETRY_VARIANT}
        voicexx::fluidsolver_${GEOMETRY_VARIANT}
        voicexx::geometry_neutrals_${GEOMETRY_VARIANT}
        voicexx::plasmasources_${GEOMETRY_VARIANT}
        voicexx::reactionrates_${GEOMETRY_VARIANT}
        voicexx::time_integration_hybrid_${GEOMETRY_VARIANT}

        paraconf::paraconf
        PDI::PDI_C
)

set_target_properties(neutrals_owngrid_${GEOMETRY_VARIANT}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
)

endforeach()
