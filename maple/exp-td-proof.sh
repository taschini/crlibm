#/bin/bash


GAPPA=~/sangria/gappaCVS/gappa/src/gappa

sed -f ./TEMPEXP/exp-td-accurate.sed exp-td-accurate1.gappa | $GAPPA 
sed -f ./TEMPEXP/exp-td-accurate.sed exp-td-accurate2.gappa | $GAPPA 
sed -f ./TEMPEXP/exp-td-accurate.sed exp-td-accurate3.gappa | $GAPPA 
sed -f ./TEMPEXP/exp-td-accurate.sed exp-td-accurate4.gappa | $GAPPA 

