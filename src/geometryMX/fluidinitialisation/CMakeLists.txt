# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("fluidinitialisation_${GEOMETRY_VARIANT}" INTERFACE)

target_include_directories("fluidinitialisation_${GEOMETRY_VARIANT}"
    INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)

target_link_libraries("fluidinitialisation_${GEOMETRY_VARIANT}" INTERFACE
        DDC::core
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::utils

)

add_library("voicexx::fluidinitialisation_${GEOMETRY_VARIANT}" ALIAS "fluidinitialisation_${GEOMETRY_VARIANT}")

endforeach()
