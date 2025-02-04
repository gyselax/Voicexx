# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("initialization_neutrals_${GEOMETRY_VARIANT}" STATIC
    restartinitializationwithneutrals.cpp
)

target_include_directories("initialization_neutrals_${GEOMETRY_VARIANT}"
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialization_neutrals_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::speciesinfo
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::utils
)

add_library("gslx::initialization_neutrals_${GEOMETRY_VARIANT}" ALIAS "initialization_neutrals_${GEOMETRY_VARIANT}")

endforeach()
