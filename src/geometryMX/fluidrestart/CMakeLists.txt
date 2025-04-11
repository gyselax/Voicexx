# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("fluidrestart_${GEOMETRY_VARIANT}" INTERFACE)

target_include_directories("fluidrestart_${GEOMETRY_VARIANT}"
    INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)

target_link_libraries("fluidrestart_${GEOMETRY_VARIANT}"
    INTERFACE
        DDC::core
        DDC::pdi
        gslx::speciesinfo
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::utils
)

add_library("voicexx::fluidrestart_${GEOMETRY_VARIANT}" ALIAS "fluidrestart_${GEOMETRY_VARIANT}")

endforeach()
