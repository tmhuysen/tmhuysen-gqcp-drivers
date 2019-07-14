# top-level CMakeLists.txt

cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

project(tmhuysen-drivers LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

message( STATUS  ${CMAKE_MODULE_PATH})
# Establish the CMake environment - define all CMake variables

find_package(gqcp 0.2.0 NO_MODULE REQUIRED)
find_package(HDF5 COMPONENTS C HL NO_MODULE REQUIRED static)

include(${CMAKE_SOURCE_DIR}/cmake/Definitions/CMakeDefinitions.cmake)
add_subdirectory(${PROJECT_DRIVERS_FOLDER})



message( STATUS  ${gqcp_INCLUDE_DIRS})