#!/bin/bash

source /usr/share/modules/init/bash

module purge
module load pgi/17.10
module load netcdf hdf5 szip zlib openmpi/pgi_2.1.2

echo "mpif90=`which mpif90`"

export FC=mpif90
export CC=mpicc
export FFLAGS="-O3 -Mvect=nosse -ta=nvidia,cc35,cuda9.0,managed,ptxinfo,fastmath,fma,unroll -Minfo=accel"
export CFLAGS="-O3"
export FREEFLAGS="-Mextend"
export FIXEDFLAGS=""
export LDFLAGS="`nf-config --flibs` `nc-config --libs` -ta=nvidia,cc35,cuda9.0,managed,ptxinfo,fastmath,fma,unroll"
export INCLUDE="-I$NETCDF_DIR/include"

export CPPDEFS=" -DFORTRANUNDERSCORE -D_SAMP_LIMIT=1000 "

