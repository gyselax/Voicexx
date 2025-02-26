# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("gridneutralfluidsolver_${GEOMETRY_VARIANT}" STATIC
    pressurediffusionsolver.cpp
    nullgridneutralfluidsolver.cpp
)

target_include_directories("gridneutralfluidsolver_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("gridneutralfluidsolver_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        sll::SLL
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::interpolation
        gslx::quadrature
        gslx::speciesinfo
        gslx::timestepper
        gslx::utils
        gslx::rhs_${GEOMETRY_VARIANT}
        voicexx::reactionrates_${GEOMETRY_VARIANT}
        voicexx::geometry_neutrals_${GEOMETRY_VARIANT}
)

add_library("voicexx::gridneutralfluidsolver_${GEOMETRY_VARIANT}" ALIAS "gridneutralfluidsolver_${GEOMETRY_VARIANT}")

endforeach()
