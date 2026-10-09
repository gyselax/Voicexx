# SPDX-License-Identifier: MIT


include(GoogleTest)

set(BASIC_GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx)

foreach(GEOMETRY_VARIANT IN LISTS BASIC_GEOMETRY_XVx_VARIANTS_LIST)

add_executable(unit_tests_reactionrates_${GEOMETRY_VARIANT}
    ../main.cpp
    reactionrates.cpp
)

target_link_libraries(unit_tests_reactionrates_${GEOMETRY_VARIANT}
    PUBLIC
        GTest::gtest
        GTest::gmock
        gyselalibxx::speciesinfo
        gyselalibxx::utils_${GEOMETRY_VARIANT}
        voicexx::fluidinitialisation_${GEOMETRY_VARIANT}
        voicexx::reactionrates_${GEOMETRY_VARIANT}
	voicexx::spline_definitions
)

gtest_discover_tests(unit_tests_reactionrates_${GEOMETRY_VARIANT}
    TEST_SUFFIX "_${GEOMETRY_VARIANT}"
    PROPERTIES TIMEOUT 10
    DISCOVERY_MODE PRE_TEST
)

endforeach()
