# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("neutralsolver_${GEOMETRY_VARIANT}" INTERFACE)

target_include_directories("neutralsolver_${GEOMETRY_VARIANT}"
    INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)

target_link_libraries("neutralsolver_${GEOMETRY_VARIANT}"
    INTERFACE
        DDC::core
        DDC::pdi
        gslx::speciesinfo
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::utils
        voicexx::fluidsolver_${GEOMETRY_VARIANT}
        voicexx::plasmaneutralscoupling_${GEOMETRY_VARIANT}
)

add_library("voicexx::neutralsolver_${GEOMETRY_VARIANT}" ALIAS "neutralsolver_${GEOMETRY_VARIANT}")

endforeach()
