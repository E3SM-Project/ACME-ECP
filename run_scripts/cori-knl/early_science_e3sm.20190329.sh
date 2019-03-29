#!/bin/bash
# Run script for SP-E3SM early science runs on summit
# 
# Branch for this simulation campaign:
# E3SM-Project/ACME-ECP.git/crjones/crm/summit_early_science (commit fdee19e)
# 
# See https://confluence.exascaleproject.org/x/SIOWAw for further details 
#
# This should serve as a template. For these summit early science 
# runs, the preferred approach is to create a new run script titled
#              early_science.$datestamp.sh 
# that corresponds to the appropriately datestamp in the case name.
#
# Contact: christopher.jones@pnnl.gov
#

create_newcase=true
dosetup=true
dobuild=true
donamelist=true
dosubmit=false

datestamp=20190329

### BASIC INFO FOR create_newcase
compset=FC5AV1C-H01A
resolution=ne120_ne120
project=m3312
machine=cori-knl
pecount=L

### Create case_name:
case_name=earlyscience.${compset}.${resolution%_*}.E3SM.$datestamp

### local directory info
repo_dir=$HOME/git_repos/Cori-Early-Science
case_dir_root=$repo_dir/Cases
case_dir=$case_dir_root/$case_name
cime_dir=$repo_dir/cime/scripts

### create case:
if [ "$create_newcase" = true ] ; then
    cd $case_dir_root
    $cime_dir/create_newcase -compset $compset -res $resolution -project $project -mach $machine -case $case_name -pecount $pecount
fi

cd $case_dir
### case setup:
if [ "$dosetup" = true ] ; then
    ./case.setup
fi

### build options:
./xmlchange ATM_PIO_NETCDF_FORMAT="64bit_data"   # note: this may not actually do anything ...
if [ "$dobuild" = true ] ; then
    ./case.build
fi

### Run options
# ./xmlchange ATM_NCPL=288
./xmlchange RUN_STARTDATE=0001-01-01
./xmlchange STOP_OPTION=nmonths
./xmlchange STOP_N=6
./xmlchange REST_OPTION=nmonths
./xmlchange REST_N=2
./xmlchange RESUBMIT=1

### namelist options
if [ "$donamelist" = true ] ; then
rm user_nl_cam
cat <<EOF >> user_nl_cam
  ! prescribed aerosol options
  prescribed_aero_cycle_yr = 01
  prescribed_aero_file = 'mam4_0.9x1.2_L72_2000clim_c170323.nc'
  prescribed_aero_datapath = '/project/projectdirs/acme/inputdata/atm/cam/chem/presc_aero'
  use_hetfrz_classnuc = .false.
  prescribed_aero_type = 'CYCLICAL'
  aerodep_flx_type = 'CYCLICAL'
  aerodep_flx_datapath = '/project/projectdirs/acme/inputdata/atm/cam/chem/presc_aero'
  aerodep_flx_file = 'mam4_0.9x1.2_L72_2000clim_c170323.nc'
  aerodep_flx_cycle_yr = 01

  ! file i/o
  nhtfrq = 0,-1,-3
  mfilt = 1,120,40
  avgflag_pertape = 'A','A','A'

  ! hourly 2D fields
  fincl2 = 'PRECT','TMQ','LHFLX','SHFLX','TS','PS','FLNT','FSNT','FSNS','FLNS','SWCF','LWCF','TGCLDLWP','TGCLDIWP'

  ! daily 3D fields + budget terms (3-hourly even better for tropical/mcs dynamics)
  fincl3 = 'T','Q','Z3','U','V','OMEGA','CLDLIQ','CLDICE','QRL','QRS'

  inithist = 'ENDOFRUN'

EOF
fi


### batch options
./xmlchange JOB_WALLCLOCK_TIME=24:00:00
./xmlchange CHARGE_ACCOUNT=$project

### submit
if [ "$dosubmit" = true ] ; then
    ./case.submit
fi
