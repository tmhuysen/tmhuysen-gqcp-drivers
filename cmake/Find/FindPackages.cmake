# Find all packages


find_package(Git REQUIRED)

find_package(Eigen3 3.3.4 REQUIRED)
find_package(gqcp 0.2.0 NO_MODULE REQUIRED)

set(BLA_VENDOR Intel10_64lp)
find_package(BLAS REQUIRED)


find_package(HDF5 COMPONENTS C HL NO_MODULE REQUIRED static)