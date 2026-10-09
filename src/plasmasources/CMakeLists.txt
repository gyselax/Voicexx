# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("plasmasources_${GEOMETRY_VARIANT}" STATIC
    maxwellian_source.cpp
)

target_include_directories("plasmasources_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("plasmasources_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gyselalibxx::geometry_${GEOMETRY_VARIANT}
        gyselalibxx::paraconfpp
        gyselalibxx::rhs_${GEOMETRY_VARIANT}
        gyselalibxx::speciesinfo
        gyselalibxx::utils

)

add_library("voicexx::plasmasources_${GEOMETRY_VARIANT}" ALIAS "plasmasources_${GEOMETRY_VARIANT}")

endforeach()
