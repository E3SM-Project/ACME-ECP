#!/bin/bash

source /sw/summitdev/lmod/7.4.0/rhel7.2_gnu4.8.5/lmod/7.4/init/sh

module purge
module load python/3.5.2
module load subversion/1.9.3
module load cmake/3.6.1
module load git/2.13.0 
module rm   xl
module load pgi/17.4
module load essl/5.5.0-20161110
module load netlib-lapack/3.6.1
module load spectrum_mpi/10.1.0.2-20161221

module load cuda/8.0.54
module load netcdf/4.4.1
module load netcdf-fortran/4.4.4
module load parallel-netcdf/1.7.0

export MPSTKZ=128M
export OMP_STACKSIZE=128M

export FC=mpif90
export CC=mpicc
export FFLAGS="-Mstack_arrays -Mextend -O2 -r8 -I${OLCF_NETCDF_FORTRAN_ROOT}/include"
export CFLAGS="-Mstack_arrays          -O2    "
export FREEFLAGS="-Mextend"
export FIXEDFLAGS=""
export LDFLAGS="-L${OLCF_NETCDF_FORTRAN_ROOT}/lib -L${OLCF_NETCDF_ROOT}/lib -lnetcdf -lnetcdff"
export INCLUDE=""

export CPPDEFS=" -DFORTRANUNDERSCORE -D_SAMP_LIMIT=512"
#export CPPDEFS="   -DSP_DIR_NS    -DCO2A -DMAXPATCH_PFT=numpft+1 -DLSMLAT=1 -DLSMLON=1 -D_MPDATA -Dm2005 -DYES3DVAL=0 -DCRM_NX=32 -DCRM_NY=1 -DCRM_NZ=58 -DCRM_DX=1000 -DCRM_DT=10  -DCRM  -DPLON=13826 -DPLAT=1 -DNP=4 -DNC=4 -DHAVE_F2003_PTR_BND_REMAP -DPLEV=72 -DPCNST=40 -DPCOLS=16 -DPSUBCOLS=1 -DN_RAD_CNST=30 -DPTRM=1 -DPTRN=1 -DPTRK=1 -D_MPI -DCAM  -D_WK_GRAD -D_PRIM  -DSPMD -DMODAL_AERO -DMODAL_AERO_4MODE_MOM   -DMODAL_AER  -DRAIN_EVAP_TO_COARSE_AERO  -DCLUBB_SGS -DCLUBB_CAM -DNO_LAPACK_ISNAN -DCLUBB_REAL_TYPE=dp -DHAVE_VPRINTF -DHAVE_TIMES -DHAVE_GETTIMEOFDAY -DHAVE_COMM_F2C -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DCNL -DCPRPGI -DNDEBUG -DMCT_INTERFACE -DHAVE_MPI -DFORTRANUNDERSCORE -DNO_SHR_VMATH -DNO_R16    -DLINUX -DCRM_DUMP -DCRM_DUMP_RATIO=0.001 -DHAVE_SLASHPROC -DUSE_CONTIGUOUS"

