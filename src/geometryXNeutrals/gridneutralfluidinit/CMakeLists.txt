# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("gridneutralfluidinit_${GEOMETRY_VARIANT}" STATIC
    gridneutralconstantinit.cpp
)

target_include_directories("gridneutralfluidinit_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("gridneutralfluidinit_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        voicexx::geometry_neutrals_${GEOMETRY_VARIANT}

)

add_library("voicexx::gridneutralfluidinit_${GEOMETRY_VARIANT}" ALIAS "gridneutralfluidinit_${GEOMETRY_VARIANT}")

endforeach()
