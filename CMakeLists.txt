
cmake_minimum_required(VERSION 3.20)
project(voicexx C CXX)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${CMAKE_SOURCE_DIR}/install")
endif()

## Look for a pre-installed paraconf
find_package(paraconf REQUIRED COMPONENTS C)

## Use the kokkos-tools for profiling purpose from `vendor/`
add_compile_definitions(KP_KERNEL_TIMER_PATH="${CMAKE_CURRENT_BINARY_DIR}/kokkos-tools/profiling/simple-kernel-timer/libkp_kernel_timer.so")


set(GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx xperiod_vx_user_defined_grid xnonperiod_vx_user_defined_grid)

set(BUILD_TESTING false)
add_subdirectory(gyselalibxx)
add_subdirectory(simulations)
