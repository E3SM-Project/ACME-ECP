#!/bin/bash

source /opt/modules/default/init/bash
module use /ccs/home/norton/.modules

module load   python/2.7.9
module rm     subversion
module load   subversion/1.8.3
module rm     cmake
module load   cmake/2.8.10.2
module rm     PrgEnv-cray
module rm     PrgEnv-gnu
module rm     PrgEnv-intel
module rm     PrgEnv-pathscale
module load   PrgEnv-pgi
module switch pgi pgi/17.1.0
module rm     cray-mpich
module rm     cray-libsci
module rm     atp
module rm     esmf
module rm     cudatoolkit
module load   cray-mpich/7.4.0
module load   cray-libsci/16.06.1
module load   atp/2.0.2
module load   esmf/5.2.0rp2
module load   cudatoolkit
module load   cray-netcdf-hdf5parallel/4.4.1.1
module load   cray-parallel-netcdf/1.7.0

export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_CPUMASK_DISPLAY=1
export MPSTKZ=128M
export OMP_STACKSIZE=128M
export CRAY_CPU_TARGET=istanbul
export CRAY_CUDA_MPS=1
export CRAY_CPU_TARGET=istanbul

export FC=ftn
export CC=cc
export FFLAGS="-O2"
export CFLAGS="-O2"
export FREEFLAGS="-Mextend"
export FIXEDFLAGS=""

