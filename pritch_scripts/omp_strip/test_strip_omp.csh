#!/bin/csh
cp control.txt test.txt
echo "BEFORE:"
cat test.txt
echo "1"
/bin/sed -i '/\!\$omp/d' ./test.txt
echo "2"
sed -i '/\!\$OMP/d' ./test.txt
echo "--------"
echo "AFTER:"
cat test.txt

