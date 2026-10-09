# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("reactionrates_${GEOMETRY_VARIANT}" STATIC
    charge_exchange.cpp
    ionisation.cpp
    recombination.cpp
    constantrate.cpp
)

target_include_directories("reactionrates_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("reactionrates_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gyselalibxx::geometry_${GEOMETRY_VARIANT}
        gyselalibxx::speciesinfo
        gyselalibxx::utils

)

add_library("voicexx::reactionrates_${GEOMETRY_VARIANT}" ALIAS "reactionrates_${GEOMETRY_VARIANT}")

endforeach()
