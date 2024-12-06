
cmake_minimum_required(VERSION 3.20)
project(voicexx C CXX)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${CMAKE_SOURCE_DIR}/install")
endif()

option(VOICEXX_BUILD_TESTING "Build the tests." ON)
option(ACTIVATE_RESTART_TESTS "Activate tests which check that a simulation gives the same results after restart." ON)

## Look for a pre-installed paraconf
find_package(paraconf REQUIRED COMPONENTS C)

add_compile_definitions(KP_KERNEL_TIMER_PATH="${CMAKE_CURRENT_BINARY_DIR}/gyselalibxx/kokkos-tools/profiling/simple-kernel-timer/libkp_kernel_timer.so")

set(GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx xperiod_vx_user_defined_grid xnonperiod_vx_user_defined_grid)

if (${VOICEXX_BUILD_TESTING})
    set(GYSELALIBXX_INCLUDE_TESTING_DEPENDENCIES ON CACHE BOOL "Ensure gyselalibxx testing dependencies are available")
endif()
set(GYSELALIBXX_BUILD_TESTING OFF)
set(GYSELALIBXX_BUILD_SIMULATIONS OFF)
add_subdirectory(gyselalibxx)

add_subdirectory(src)
add_subdirectory(simulations)


if (${VOICEXX_BUILD_TESTING})
    find_package(GTest "1.12" QUIET)
    ## Use CTest for running tests
    include(CTest)
    add_subdirectory(tests)
endif()
