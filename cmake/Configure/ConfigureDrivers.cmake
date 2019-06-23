# Configure and install the drivers

message(STATUS ${PROJECT_DRIVERS_SOURCE_FILES})


foreach(DRIVER_SOURCE ${PROJECT_DRIVERS_SOURCE_FILES})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(DRIVER_NAME ${DRIVER_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${DRIVER_NAME} ${DRIVER_SOURCE})

    # Include gqcp
    target_include_directories(${DRIVER_NAME} PUBLIC ${gqcp_INCLUDE_DIRS})
    target_link_libraries(${DRIVER_NAME} PUBLIC gqcp)
endforeach()
