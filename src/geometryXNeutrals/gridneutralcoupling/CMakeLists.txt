# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("gridneutralcoupling_${GEOMETRY_VARIANT}" STATIC
    nullgridneutralcoupling.cpp
)

target_include_directories("gridneutralcoupling_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("gridneutralcoupling_${GEOMETRY_VARIANT}"
    PUBLIC
        gslx::geometry_${GEOMETRY_VARIANT}
        voicexx::geometry_neutrals_${GEOMETRY_VARIANT}
)

add_library("voicexx::gridneutralcoupling_${GEOMETRY_VARIANT}" ALIAS "gridneutralcoupling_${GEOMETRY_VARIANT}")

endforeach()
