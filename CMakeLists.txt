# top-level CMakeLists.txt

cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

project(tmhuysen-drivers LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

# Establish the CMake environment - define all CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/Definitions/CMakeDefinitions.cmake)
add_subdirectory(${PROJECT_DRIVERS_FOLDER})

message( STATUS  ${gqcp_INCLUDE_DIRS})