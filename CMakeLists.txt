# top-level CMakeLists.txt

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(tmhuysen-drivers LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

# Establish the CMake environment - define all CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/Definitions/CMakeDefinitions.cmake)
add_subdirectory(${PROJECT_DRIVERS_FOLDER})
