# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("gridneutraltimeintegration_${GEOMETRY_VARIANT}" STATIC
    gridneutralpredcorr.cpp
)

target_include_directories("gridneutraltimeintegration_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("gridneutraltimeintegration_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::initialisation_${GEOMETRY_VARIANT}
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::utils
        voicexx::geometry_neutrals_${GEOMETRY_VARIANT}
        voicexx::plasmaneutralscoupling_${GEOMETRY_VARIANT}
        voicexx::gridneutralfluidsolver_${GEOMETRY_VARIANT}

)

add_library("voicexx::gridneutraltimeintegration_${GEOMETRY_VARIANT}" ALIAS "gridneutraltimeintegration_${GEOMETRY_VARIANT}")

endforeach()
