# SPDX-License-Identifier: MIT


foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("geometry_neutrals_${GEOMETRY_VARIANT}" STATIC
    gridneutral_interpolator.cpp
)
target_include_directories("geometry_neutrals_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries("geometry_neutrals_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        gyselalibxx::speciesinfo
        gyselalibxx::moments
        gyselalibxx::utils
        gyselalibxx::geometry_${GEOMETRY_VARIANT}
	voicexx::spline_definitions
)
add_library("voicexx::geometry_neutrals_${GEOMETRY_VARIANT}" ALIAS "geometry_neutrals_${GEOMETRY_VARIANT}")

endforeach()
