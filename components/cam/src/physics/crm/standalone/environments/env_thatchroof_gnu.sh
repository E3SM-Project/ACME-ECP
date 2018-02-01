#!/bin/bash

source /usr/share/modules/init/bash

module purge
module load netcdf hdf5 szip zlib curl

echo "mpif90=`which mpif90`"

export FC=mpif90
export CC=mpicc
export FFLAGS="-O2"
export CFLAGS="-O2"
export FREEFLAGS="-ffree-line-length-none -ffixed-line-length-none"
export FIXEDFLAGS=""
export LDFLAGS="`nf-config --flibs` `nc-config --libs`"
export INCLUDE="-I$NETCDF_DIR/include"

export CPPDEFS=" -DFORTRANUNDERSCORE "

