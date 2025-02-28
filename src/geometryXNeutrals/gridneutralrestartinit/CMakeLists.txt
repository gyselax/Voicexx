# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("gridneutralrestartinit_${GEOMETRY_VARIANT}" STATIC
    gridneutralrestartinit.cpp
)

target_include_directories("gridneutralrestartinit_${GEOMETRY_VARIANT}"
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("gridneutralrestartinit_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::utils
        voicexx::geometry_neutrals_${GEOMETRY_VARIANT}
)

add_library("voicexx::gridneutralrestartinit_${GEOMETRY_VARIANT}" ALIAS "gridneutralrestartinit_${GEOMETRY_VARIANT}")

endforeach()
