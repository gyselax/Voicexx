
cmake_minimum_required(VERSION 3.20)
project(voicexx C CXX)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${CMAKE_SOURCE_DIR}/install")
endif()

option(VOICEXX_BUILD_TESTING "Build the tests." ON)

## Look for a pre-installed paraconf
find_package(paraconf REQUIRED COMPONENTS C)


set(GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx xperiod_vx_user_defined_grid xnonperiod_vx_user_defined_grid)

set(GYSELALIBXX_INCLUDE_TESTING_DEPENDENCIES ON)
set(GYSELALIBXX_BUILD_TESTING OFF)
set(GYSELALIBXX_BUILD_SIMULATIONS OFF)
add_subdirectory(gyselalibxx)

add_subdirectory(src)
add_subdirectory(simulations)


if (${VOICEXX_BUILD_TESTING})
add_subdirectory(tests)
endif()
