# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("initialisation_neutrals_${GEOMETRY_VARIANT}" STATIC
    restartinitialisationwithneutrals.cpp
)

target_include_directories("initialisation_neutrals_${GEOMETRY_VARIANT}"
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialisation_neutrals_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::speciesinfo
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::utils
)

add_library("voicexx::initialisation_neutrals_${GEOMETRY_VARIANT}" ALIAS "initialisation_neutrals_${GEOMETRY_VARIANT}")

endforeach()
