# Define the CMake environment

# Uppercase and lowercase names
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)

# Specify the library build type
set(LIBRARY_NAME ${PROJECT_NAME})
set(LIBRARY_TYPE SHARED)  # dynamic/shared library

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Building ${LIBRARY_NAME} in ${CMAKE_BUILD_TYPE} mode")

include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Compiler.cmake)

include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Options.cmake)

# Include driver and source files
include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Drivers.cmake)

# Find all packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Find)  # let CMake know that we have supplied our own FindXXX.cmake files
include(${CMAKE_SOURCE_DIR}/cmake/Find/FindPackages.cmake)

include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Git.cmake)
