# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("reactionrates_${GEOMETRY_VARIANT}" STATIC
    constantrate.cpp
    charge_exchange.cpp
    ionization.cpp
    recombination.cpp
)

target_include_directories("reactionrates_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("reactionrates_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        sll::SLL
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::utils

)

add_library("voicexx::reactionrates_${GEOMETRY_VARIANT}" ALIAS "reactionrates_${GEOMETRY_VARIANT}")

endforeach()
