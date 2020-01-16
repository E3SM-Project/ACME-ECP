#!/bin/csh
set sourceroot = "/home1/00993/tg802402/repositories/ACME-ECP/components"
foreach file ( ` find $sourceroot -name "*.[fF]*" ` )
  echo $file
end
