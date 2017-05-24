#!/bin/bash

export FC=mpif90
export CC=mpicc
export FFLAGS="-O2"
export CFLAGS="-O2"
export FREEFLAGS="-ffree-line-length-none"
export FIXEDFLAGS=""
export LDFLAGS="${NCLIBS}"
export INCLUDE="${NCINC}"

export CPPDEFS=" -DFORTRANUNDERSCORE "

