# List all drivers


# Find the drivers folder
set(PROJECT_DRIVERS_FOLDER ${CMAKE_SOURCE_DIR}/drivers)


# Find the source files for the drivers
set(PROJECT_DRIVERS_SOURCE_FILES
        ${PROJECT_DRIVERS_FOLDER}/davidson_constrained_N_O.cpp
        ${PROJECT_DRIVERS_FOLDER}/dense_constrained_N_O.cpp
        ${PROJECT_DRIVERS_FOLDER}/test_driver.cpp)
