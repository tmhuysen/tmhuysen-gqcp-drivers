# Configure and install the drivers

message(STATUS ${PROJECT_DRIVERS_SOURCE_FILES})


foreach(DRIVER_SOURCE ${PROJECT_DRIVERS_SOURCE_FILES})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(DRIVER_NAME ${DRIVER_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${DRIVER_NAME} ${DRIVER_SOURCE})

    target_include_directories(${DRIVER_NAME}
            PUBLIC
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
            PRIVATE
            ${CMAKE_SOURCE_DIR}/src
            )

    target_link_libraries(${DRIVER_NAME}
            PUBLIC
            Eigen3::Eigen
            gqcp::gqcp
            )


    target_include_directories(${DRIVER_NAME} PUBLIC $<BUILD_INTERFACE:${BLAS_INCLUDE_DIR}>)
    target_link_libraries(${DRIVER_NAME} PUBLIC ${BLAS_LIBRARIES})
    target_compile_options(${DRIVER_NAME} PUBLIC -DEIGEN_USE_MKL_ALL -DMKL_LP64)
endforeach()
