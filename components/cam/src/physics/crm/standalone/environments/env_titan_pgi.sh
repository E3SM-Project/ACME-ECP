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
module switch pgi pgi/17.3.0
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
module load   ncl
module load   nco

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
export FFLAGS="-Mstack_arrays -Mextend -O2 -r8"
export CFLAGS="-Mstack_arrays          -O2    "
export FREEFLAGS="-Mextend"
export FIXEDFLAGS=""
export LDFLAGS=""
export INCLUDE=""

export CPPDEFS=" -DFORTRANUNDERSCORE -D_SAMP_LIMIT=512"
#export CPPDEFS="   -DSP_DIR_NS    -DCO2A -DMAXPATCH_PFT=numpft+1 -DLSMLAT=1 -DLSMLON=1 -D_MPDATA -Dm2005 -DYES3DVAL=0 -DCRM_NX=32 -DCRM_NY=1 -DCRM_NZ=58 -DCRM_DX=1000 -DCRM_DT=10  -DCRM  -DPLON=13826 -DPLAT=1 -DNP=4 -DNC=4 -DHAVE_F2003_PTR_BND_REMAP -DPLEV=72 -DPCNST=40 -DPCOLS=16 -DPSUBCOLS=1 -DN_RAD_CNST=30 -DPTRM=1 -DPTRN=1 -DPTRK=1 -D_MPI -DCAM  -D_WK_GRAD -D_PRIM  -DSPMD -DMODAL_AERO -DMODAL_AERO_4MODE_MOM   -DMODAL_AER  -DRAIN_EVAP_TO_COARSE_AERO  -DCLUBB_SGS -DCLUBB_CAM -DNO_LAPACK_ISNAN -DCLUBB_REAL_TYPE=dp -DHAVE_VPRINTF -DHAVE_TIMES -DHAVE_GETTIMEOFDAY -DHAVE_COMM_F2C -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DCNL -DCPRPGI -DNDEBUG -DMCT_INTERFACE -DHAVE_MPI -DFORTRANUNDERSCORE -DNO_SHR_VMATH -DNO_R16    -DLINUX -DCRM_DUMP -DCRM_DUMP_RATIO=0.001 -DHAVE_SLASHPROC -DUSE_CONTIGUOUS"

