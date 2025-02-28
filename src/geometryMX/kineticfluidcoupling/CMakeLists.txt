# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("fluidrhs_${GEOMETRY_VARIANT}" STATIC
    kinetic_fluid_coupling_source.cpp
)

target_include_directories("fluidrhs_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("fluidrhs_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::initialisation_${GEOMETRY_VARIANT}
        gslx::quadrature
        gslx::rhs_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::timestepper
        gslx::utils
        gslx::utils_${GEOMETRY_VARIANT}
        voicexx::reactionrates_${GEOMETRY_VARIANT}
)

add_library("voicexx::fluidrhs_${GEOMETRY_VARIANT}" ALIAS "fluidrhs_${GEOMETRY_VARIANT}")

endforeach()
