# CMake initial cache file for titan
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF_DIR} CACHE FILEPATH "")
SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (HDF5_DIR $ENV{HDF5_DIR} CACHE FILEPATH "")
SET (SZIP_DIR $ENV{SZIP_DIR} CACHE FILEPATH "")


SET (USE_QUEUING FALSE CACHE BOOL "")