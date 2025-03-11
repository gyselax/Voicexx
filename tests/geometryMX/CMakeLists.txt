# SPDX-License-Identifier: MIT

find_package(Python3 REQUIRED COMPONENTS Interpreter)

include(GoogleTest)

set(BASIC_GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx)

foreach(GEOMETRY_VARIANT IN LISTS BASIC_GEOMETRY_XVx_VARIANTS_LIST)

add_executable(unit_tests_geometryMX_${GEOMETRY_VARIANT}
	../main.cpp
    fluidspecies.cpp
    moments.cpp
    predcorr_hybrid.cpp
    diffusiveneutrals.cpp
    kineticfluidcoupling.cpp
)

target_link_libraries(unit_tests_geometryMX_${GEOMETRY_VARIANT}
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::advection
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::moments
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::quadrature
        gslx::speciesinfo
        gslx::time_integration_${GEOMETRY_VARIANT}
        gslx::utils_${GEOMETRY_VARIANT}
        voicexx::fluidinitialisation_${GEOMETRY_VARIANT}
        voicexx::fluidsolver_${GEOMETRY_VARIANT}
        voicexx::time_integration_hybrid_${GEOMETRY_VARIANT}
        voicexx::reactionrates_${GEOMETRY_VARIANT}
)

gtest_discover_tests(unit_tests_geometryMX_${GEOMETRY_VARIANT}
    TEST_SUFFIX "_${GEOMETRY_VARIANT}"
    PROPERTIES TIMEOUT 10
    DISCOVERY_MODE PRE_TEST
)

endforeach()

if (${ACTIVATE_RESTART_TESTS})
    add_test(NAME TestSimulationNeutralsRestart_xnonperiod_vx
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_neutrals_restart.sh"
            "${PROJECT_SOURCE_DIR}"
            "$<TARGET_FILE:neutrals_xnonperiod_vx>"
            "$<TARGET_FILE:Python3::Interpreter>"
            "restart"
            "1e-15"
            "1e-15")
        set_property(TEST TestSimulationNeutralsRestart_xnonperiod_vx PROPERTY TIMEOUT 200)
        set_property(TEST TestSimulationNeutralsRestart_xnonperiod_vx PROPERTY COST 100)
endif()
