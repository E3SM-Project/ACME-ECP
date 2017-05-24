#!/bin/bash

function usage {
  echo "./compute_diffs.sh -f1 </path/file1> -f2 </path/file2>"
  echo "REQUIRED:"
  echo "   -f1 </path/file1>             First file"
  echo "   -f2 </path/file2>             Second file"
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
    -f1)
    F1="$2"
    shift # past argument
    ;;
    -f2)
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

[ "$F1" == "" ] && echo "Must specify -f1" && exit -1
[ "$F2" == "" ] && echo "Must specify -f2" && exit -1

ncl -h || echo "NCL not installed / loaded."
ncl -h || exit -1

ncl "fn1=\"$F1\"" "fn2=\"$F2\"" diff.ncl


