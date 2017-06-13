#!/bin/bash

#Just make sure you have NetCDF, NCO, and NCL installed, and you define NCLIBS and NCINC
#with all needed link line libraries and includes, respectively

export FC=mpif90
export CC=mpicc
export FFLAGS="-O2"
export CFLAGS="-O2"
export FREEFLAGS="-ffree-line-length-none"
export FIXEDFLAGS=""
export LDFLAGS="${NCLIBS}"
export INCLUDE="${NCINC}"

export CPPDEFS=" -DFORTRANUNDERSCORE -D_SAMP_LIMIT=32 "

