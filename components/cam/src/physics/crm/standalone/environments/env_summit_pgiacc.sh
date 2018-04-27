#!/bin/bash

source $MODULESHOME/init/bash

module purge
module load pgi/18.3
module load spectrum-mpi
module load netcdf netcdf-fortran
module load cuda

echo "mpif90=`which mpif90`"

export FC=mpif90
export CC=mpicc
export FFLAGS="-O3 -Mvect=nosimd -ta=nvidia,cc70,cuda9.0,pinned,ptxinfo -Minfo=accel"
export CFLAGS="-O3"
export FREEFLAGS="-Mextend"
export FIXEDFLAGS=""
export LDFLAGS="`nf-config --flibs` `nc-config --libs` -ta=nvidia,cc70,cuda9.0,pinned,ptxinfo"
export INCLUDE="-I$OLCF_NETCDF_ROOT/include -I$OLCF_NETCDF_FORTRAN_ROOT/include"

export CPPDEFS=" -DFORTRANUNDERSCORE "

