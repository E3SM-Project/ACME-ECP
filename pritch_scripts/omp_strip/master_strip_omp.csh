#!/bin/csh
set sourceroot = "/home1/00993/tg802402/repositories/ACME-ECP/components"
foreach file ( ` find $sourceroot -name "*.[fF]*" ` )
sed -i '/\!\$omp/d' $file
sed -i '/\!\$OMP/d' $file
echo "Done $file"
end

