#!/bin/bash

function usage {
  echo "./compute_diffs.sh -p1 </path/file1> -p2 </path/file2>"
  echo "REQUIRED:"
  echo "   -p1 </path/prefix1>             First file"
  echo "   -p2 </path/prefix2>             Second file"
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
    -p1)
    F1="$2"
    shift # past argument
    ;;
    -p2)
    F2="$2"
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

[ "$F1" == "" ] && echo "Must specify -p1" && usage && exit -1
[ "$F2" == "" ] && echo "Must specify -p2" && usage && exit -1

ncl -h >& /dev/null || echo "NCL not installed / loaded."
ncl -h >& /dev/null || exit -1

for f in `ls ${F1}_*.nc`; do
  suffix=${f#$F1}
  ncl "fn1=\"${F1}${suffix}\"" "fn2=\"${F2}${suffix}\"" diff.ncl
done


