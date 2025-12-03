# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("time_integration_hybrid_${GEOMETRY_VARIANT}" INTERFACE)

target_include_directories("time_integration_hybrid_${GEOMETRY_VARIANT}"
    INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)

target_link_libraries("time_integration_hybrid_${GEOMETRY_VARIANT}"
    INTERFACE
        DDC::core
        DDC::pdi
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::initialisation_${GEOMETRY_VARIANT}
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::utils
        voicexx::fluidsolver_${GEOMETRY_VARIANT}
        voicexx::neutralsolver_${GEOMETRY_VARIANT}
        voicexx::plasmaneutralscoupling_${GEOMETRY_VARIANT}

)

add_library("voicexx::time_integration_hybrid_${GEOMETRY_VARIANT}" ALIAS "time_integration_hybrid_${GEOMETRY_VARIANT}")

endforeach()
