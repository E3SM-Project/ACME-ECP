#!/bin/bash

source /usr/share/modules/init/bash

module purge
module load pgi/16.10
module load netcdf hdf5 szip zlib curl openmpi

echo "mpif90=`which mpif90`"

export FC=mpif90
export CC=mpicc
export FFLAGS="-O2 -ta=tesla,unroll,cuda8.0,cc60,loadcache:L1 -Minfo=accel"
export CFLAGS="-O2 -ta=tesla,unroll,cuda8.0,cc60,loadcache:L1 -Minfo=accel"
export FREEFLAGS="-Mextend"
export FIXEDFLAGS=""
export LDFLAGS="-L$NETCDF_DIR/lib -lnetcdff -lnetcdf -L$HDF5_DIR/lib -lhdf5_hl -lhdf5 -L$SZIP_DIR/lib -lsz -L$ZLIB_DIR/lib -lz -L$CURL_DIR/lib -lcurl -ldl -ta=tesla,pinned,unroll,cuda8.0,cc60,loadcache:L1"
export INCLUDE="-I$NETCDF_DIR/include"

export CPPDEFS=" -DFORTRANUNDERSCORE -D_SAMP_LIMIT=32 "

