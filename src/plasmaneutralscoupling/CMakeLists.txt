# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("plasmaneutralscoupling_${GEOMETRY_VARIANT}" STATIC
    noenergytransfercoupling.cpp
    densitycoupling.cpp
)

target_include_directories("plasmaneutralscoupling_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("plasmaneutralscoupling_${GEOMETRY_VARIANT}"
    PUBLIC
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::rhs_${GEOMETRY_VARIANT}
        voicexx::geometry_neutrals_${GEOMETRY_VARIANT}
        voicexx::reactionrates_${GEOMETRY_VARIANT}
)

add_library("voicexx::plasmaneutralscoupling_${GEOMETRY_VARIANT}" ALIAS "plasmaneutralscoupling_${GEOMETRY_VARIANT}")

endforeach()
