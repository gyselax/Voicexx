# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("fluidsolver_${GEOMETRY_VARIANT}" STATIC
    samegridfluidsolver.cpp
    diffgridsfluidsolver.cpp
)

target_include_directories("fluidsolver_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("fluidsolver_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::interpolation
        gslx::math_tools
        gslx::quadrature
        gslx::rhs_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::timestepper
        gslx::utils
        voicexx::geometry_neutrals_${GEOMETRY_VARIANT}
        voicexx::reactionrates_${GEOMETRY_VARIANT}
)

add_library("voicexx::fluidsolver_${GEOMETRY_VARIANT}" ALIAS "fluidsolver_${GEOMETRY_VARIANT}")

endforeach()
