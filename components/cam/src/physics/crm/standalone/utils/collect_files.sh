#!/bin/bash

function usage {
  echo "./collect_files.sh -prefix </path/file_prefix>"
  echo "REQUIRED:"
  echo "   -prefix <path>                Prefix for files to collect into a single file (ncrcat)"
  echo "   -h | --help                   Print help"
  exit 0
}

##########################################################################################
## READ IN COMMAND LINE PARAMETERS
##########################################################################################

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -prefix)
    MYPREFIX="$2"
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

[ "$MYPREFIX" == "" ] && echo "Must specify -prefix" && exit -1

ncrcat --help || echo "NCO not installed / loaded. ncrcat not available"
ncrcat --help || exit -1

ncrcat ${MYPREFIX}_*.nc ${MYPREFIX}.nc && rm ${MYPREFIX}_*.nc || exit -1


