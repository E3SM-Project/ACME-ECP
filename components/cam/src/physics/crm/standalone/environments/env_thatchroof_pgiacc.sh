#!/bin/bash

source /usr/share/modules/init/bash

module purge
module load pgi/17.10
module load netcdf hdf5 szip zlib openmpi/pgi_2.1.2

echo "mpif90=`which mpif90`"

export FC=mpif90
export CC=mpicc
export FFLAGS="-I$NETCDF_DIR/include -O3 -ta=nvidia,cc50,managed,ptxinfo -Minfo=accel"
export CFLAGS="-O3"
export FREEFLAGS="-Mextend"
export FIXEDFLAGS=""
export LDFLAGS="`nf-config --flibs` `nc-config --libs` -ta=nvidia,cc50,managed,ptxinfo"
export INCLUDE="-I$NETCDF_DIR/include"

export CPPDEFS=" -DFORTRANUNDERSCORE "

