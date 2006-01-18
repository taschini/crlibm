#/bin/bash


GAPPA=~/schorle/gappa/src/gappa

echo Proving error bound on sqrt12...
sed -f sqrt.sed sqrt12.gappa | $GAPPA
echo Proving error bound on sqrt13...
sed -f sqrt.sed sqrt13.gappa | $GAPPA

