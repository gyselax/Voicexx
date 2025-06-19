# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("geometry_moments_${GEOMETRY_VARIANT}" INTERFACE)
target_include_directories("geometry_moments_${GEOMETRY_VARIANT}"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_moments_${GEOMETRY_VARIANT}" INTERFACE
    DDC::core
    gslx::geometry_${GEOMETRY_VARIANT}
)
add_library("voicexx::geometry_moments_${GEOMETRY_VARIANT}" ALIAS "geometry_moments_${GEOMETRY_VARIANT}")

endforeach()
