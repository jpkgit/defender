cmake_minimum_required(VERSION 2.8.12)
project(defender C)
set(PACKAGE defender)
include(${PROJECT_SOURCE_DIR}/../cmake/set_release.cmake)
add_definitions(-DTOOL_RELEASE="${RELEASE}")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/modules)

if(MSVC)
include_directories(getopt)
add_definitions(/D _CRT_SECURE_NO_WARNINGS)
else()
add_definitions(-Wall)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu90")
endif()

if(NOT libhackrf_SOURCE_DIR)
find_package(LIBHACKRF REQUIRED)
include_directories(${LIBHACKRF_INCLUDE_DIR})
else()
include_directories(${libhackrf_SOURCE_DIR}/src)
endif()

add_subdirectory(src)

########################################################################
# Create uninstall target
########################################################################

if(NOT HackRF_SOURCE_DIR)
configure_file(
    ${PROJECT_SOURCE_DIR}/../cmake/cmake_uninstall.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
@ONLY)

add_custom_target(uninstall
    ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)
endif()
