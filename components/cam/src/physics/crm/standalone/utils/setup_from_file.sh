#!/bin/bash

function usage {
  echo "./setup.sh -crm_root <path> -build_root <path> -file <file>"
  echo "REQUIRED:"
  echo "   -crm_root <path>                Absolute path to the crm source root (/[...]/crm)"
  echo "   -build_root <path>              Absolute path to the build directory root"
  echo "   -file </path/to/dump/file>      Standolone CRM dump file"
  echo "   -h | --help                     Print help"
  exit 0
}

##########################################################################################
## READ IN COMMAND LINE PARAMETERS
##########################################################################################
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
    -file)
    TFILE="$2"
    shift # past argument
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
## Make sure they specified all variables
##########################################################################################
[ "$CRM_ROOT"   == "" ] && echo "Must specify -crm_root <path>" && exit -1
[ "$BUILD_ROOT" == "" ] && echo "Must specify -build_root <path>" && exit -1
[ "$TFILE"      == "" ] && echo "Must specify -file <path>" && exit -1

##########################################################################################
## Make sure NetCDF is in the PATH
##########################################################################################
ncdump >& /dev/null || echo "NetCDF is currently not in the PATH. Please source the correct environment before running this script"
ncdump >& /dev/null || exit -1

##########################################################################################
## Gather the information we want from the netcdf dump file
##########################################################################################
CRM_NX=`  ncdump -h $TFILE | grep crm_nx    | tail -n1 | cut -d ' ' -f3`
CRM_NY=`  ncdump -h $TFILE | grep crm_ny    | tail -n1 | cut -d ' ' -f3`
CRM_NZ=`  ncdump -h $TFILE | grep crm_nz    | tail -n1 | cut -d ' ' -f3`
CRM_DX=`  ncdump -h $TFILE | grep crm_dx    | tail -n1 | cut -d ' ' -f3`
CRM_DT=`  ncdump -h $TFILE | grep crm_dt    | tail -n1 | cut -d ' ' -f3`
MICRO=`   ncdump -h $TFILE | grep micro     | tail -n1 | cut -d ' ' -f3 | cut -d '"' -f 2`
ADV=`     ncdump -h $TFILE | grep crm_adv   | tail -n1 | cut -d ' ' -f3 | cut -d '"' -f 2`
PLEV=`    ncdump -h $TFILE | grep plev      | tail -n1 | cut -d ' ' -f3`
PSUBCOLS=`ncdump -h $TFILE | grep psubcols  | tail -n1 | cut -d ' ' -f3`
PCOLS=`   ncdump -h $TFILE | grep pcols     | tail -n1 | cut -d ' ' -f3`
PCNST=`   ncdump -h $TFILE | grep pcnst     | tail -n1 | cut -d ' ' -f3`

echo "./setup.sh -crm_root $CRM_ROOT -build_root $BUILD_ROOT -nadv $PCNST -nlev $PLEV -crm_nx $CRM_NX -crm_ny $CRM_NY -crm_nz $CRM_NZ -crm_dx $CRM_DX -crm_dt $CRM_DT -pcols $PCOLS -psubcols $PSUBCOLS -SPCAM_microp_scheme $MICRO -crm_adv $ADV"

./setup.sh -crm_root $CRM_ROOT -build_root $BUILD_ROOT -nadv $PCNST -nlev $PLEV -crm_nx $CRM_NX -crm_ny $CRM_NY -crm_nz $CRM_NZ -crm_dx $CRM_DX -crm_dt $CRM_DT -pcols $PCOLS -psubcols $PSUBCOLS -SPCAM_microp_scheme $MICRO -crm_adv $ADV




