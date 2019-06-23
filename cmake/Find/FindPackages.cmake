# Find all packages


find_package(Git REQUIRED)

find_package(Eigen3 3.3.4 REQUIRED)
find_package(gqcp 0.2.0 REQUIRED)

if (USE_MKL)
    find_package(MKL REQUIRED)
endif()

find_package(HDF5 COMPONENTS C HL NO_MODULE REQUIRED static)