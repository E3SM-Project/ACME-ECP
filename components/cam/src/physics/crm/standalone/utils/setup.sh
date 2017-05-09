#!/bin/bash

function usage {
  echo "./setup.sh -crm_root <path> [ ... ]"
  echo "REQUIRED:"
  echo "   -crm_root <path>                Absolute path to the crm source root (/[...]/crm)"
  echo "   -build_root <path>              Absolute path to the build directory root"
  echo "OPTIONAL:"
  echo "   -nadv <n>                       Number of advective constituents (PCNST)"
  echo "   -nlev <n>                       Number of CAM vertical levels"
  echo "   -crm_nx <n>                     CRM's x-grid."
  echo "   -crm_ny <n>                     CRM's y-grid."
  echo "   -crm_nz <n>                     CRM's z-grid."
  echo "   -crm_dx <n>                     CRM's horizontal grid spacing."
  echo "   -crm_dt <n>                     CRM's timestep."
  echo "   -pcols <n>                      PCOLS"
  echo "   -psubcols <n>                   PSUBCOLS"
  echo "   -SPCAM_microp_scheme <string>   CRM microphysics package name [sam1mom | m2005 ]."
  echo "   -clubb_crm                      CRM with clubb treatment"
  echo "   -use_crm_cldfrac                Use fractional cloudiness in CRM"
  echo "   -crm_adv                        CRM advection scheme [MPDATA | UM5]"
  echo "   -h | --help                     Print help"
  exit 0
}

##########################################################################################
## READ IN COMMAND LINE PARAMETERS
##########################################################################################

#Defaults
BUILD_ROOT=""
CRM_ROOT=""
PLEV=30
PCOLS=16
PSUBCOLS=1
PCNST=50
CRM_NX=32
CRM_NY=1
CRM_NZ=28
CRM_DX=1000
CRM_DT=5
MICRO=m2005
CLUBB_CRM=0
USE_CRM_CLDFRAC=0
ADV=MPDATA
yes3Dval=0

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -crm_root)
    CRM_ROOT="$2"
    shift # past argument
    ;;
    -build_root)
    BUILD_ROOT="$2"
    shift # past argument
    ;;
    -nadv)
    PCNST="$2"
    shift # past argument
    ;;
    -pcols)
    PCOLS="$2"
    shift # past argument
    ;;
    -psubcols)
    PSUBCOLS="$2"
    shift # past argument
    ;;
    -nlev)
    PLEV="$2"
    shift # past argument
    ;;
    -crm_nx)
    CRM_NX="$2"
    shift # past argument
    ;;
    -crm_ny)
    CRM_NY="$2"
    shift # past argument
    ;;
    -crm_nz)
    CRM_NZ="$2"
    shift # past argument
    ;;
    -crm_dx)
    CRM_DX="$2"
    shift # past argument
    ;;
    -crm_dt)
    CRM_DT="$2"
    shift # past argument
    ;;
    -SPCAM_microp_scheme)
    MICRO="$2"
    shift # past argument
    ;;
    -crm_adv)
    ADV=$2
    shift # past argument
    ;;
    -clubb_crm)
    CLUBB_CRM=1
    ;;
    -use_crm_cldfrac)
    USE_CRM_CLDFRAC=1
    ;;
    -h|--help)
    usage
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

##########################################################################################
## CHECK SOME OF THE PARAMETERS
##########################################################################################
[ "$CRM_ROOT"   == "" ] && echo "Must specify -crm_root <path> as a parameter" && exit -1
[ "$BUILD_ROOT" == "" ] && echo "Must specify -build_root <path> as a parameter" && exit -1
[ "$MICRO" != "sam1mom" ] && [ "$MICRO" != "m2005" ] && echo "-SPCAM_microp_scheme must be [m2005 | sam1mom]" && exit -1
[ "$ADV"   != "MPDATA"  ] && [ "$ADV"   != "UM5"   ] && echo "-crm_adv must be [MPDATA | UM5]" && exit -1
[ "$CRM_NY" -gt "1" ] && yes3Dval=1

mkdir -p $BUILD_ROOT
SCRIPTS=`pwd`
cd $BUILD_ROOT

##########################################################################################
## Create the preprocessor define flags
##########################################################################################
CPPDEFS=" -DCRM -D$MICRO -DYES3DVAL=$yes3Dval -DCRM_NX=$CRM_NX -DCRM_NY=$CRM_NY -DCRM_NZ=$CRM_NZ -DCRM_DX=$CRM_DX -DCRM_DT=$CRM_DT -DPLEV=$PLEV -DPSUBCOLS=$PSUBCOLS -DPCOLS=$PCOLS -DPCNST=$PCNST  -DHAVE_IEEE_ARITHMETIC -DCRM_STANDALONE"
[ "$CLUBB_CRM" -eq "1" ] && CPPDEFS="${CPPDEFS} -DCLUBB_CRM -DCLUBB_REAL_TYPE=dp "

##########################################################################################
## Determine the directories to inluce for source files
##########################################################################################
echo "$CRM_ROOT"                                               > Filepath
echo "$CRM_ROOT/standalone/src"                               >> Filepath
[ "$CLUBB_CRM" -eq "1"  ] && echo "$CRM_ROOT/../clubb"        >> Filepath \
                          && echo "$CRM_ROOT/SGS_CLUBBkvhkvm" >> Filepath
[ "$CLUBB_CRM" -ne "1"  ] && echo "$CRM_ROOT/SGS_TKE"         >> Filepath
[ "$ADV"   == "UM5"     ] && echo "$CRM_ROOT/ADV_UM5"         >> Filepath
[ "$ADV"   == "MPDATA"  ] && echo "$CRM_ROOT/ADV_MPDATA"      >> Filepath
[ "$MICRO" == "sam1mom" ] && echo "$CRM_ROOT/MICRO_SAM1MOM"   >> Filepath
[ "$MICRO" == "m2005"   ] && [ "$USE_CRM_CLDFRAC" ==  "1" ] && echo "$CRM_ROOT/MICRO_M2005FRAC" >> Filepath
[ "$MICRO" == "m2005"   ] && [ "$USE_CRM_CLDFRAC" -ne "1" ] && echo "$CRM_ROOT/MICRO_M2005"     >> Filepath

##########################################################################################
## Copy all potentially needed source files to build_root
##########################################################################################
for i in `cat Filepath`
  do cp -f $i/*.f $i/*.F $i/*.f90 $i/*.F90 $i/*.c $i/*.h $i/*.in $i/*.inc . >& /dev/null
done
cp $CRM_ROOT/../cam/crmdims.F90 .
cp $CRM_ROOT/../cam/ppgrid.F90 .
cp $CRM_ROOT/../../../../../cime/share/csm_share/shr/shr_const_mod.F90 .
cp $CRM_ROOT/../../../../../cime/share/csm_share/shr/shr_kind_mod.F90 .

##########################################################################################
## Define the dependencies between the source files
##########################################################################################
ls *.f *.F *.F90 *.f90 *.c > Srcfiles 2> /dev/null
echo "." > Filepath
$SCRIPTS/mkDepends ./Filepath ./Srcfiles > ./Depends

##########################################################################################
## WRite preprocessor directives, and copy Makefile
##########################################################################################
echo "CPPDEFS := ${CPPDEFS}" > make.inc
cp $SCRIPTS/Makefile .

