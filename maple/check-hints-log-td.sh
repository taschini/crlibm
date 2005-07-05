#!/bin/sh
#
# This will only function if you have a modified version of the gappa tool
#
# You probably need to edit the path to the gappa executable
GAPPA=~/ble/gappa-0.4.5/src/gappa


cat log-td.gappa | grep "#MAPLE" | sed -e s/=/:=/ -e "s/;/:/" > log-td.gappa.autocheck.mpl
echo "printf(\"You should see only zeros in what follows:\\n\");" >> log-td.gappa.autocheck.mpl
sed -f ./TEMPLOG/log-td_1.sed log-td.gappa | $GAPPA 2>&1 | grep simplify >> log-td.gappa.autocheck.mpl
echo "printf(\"If you have seen only zeros up to now, everything's fine \\n\");" >> log-td.gappa.autocheck.mpl

cat log-td-E0.gappa | grep "#MAPLE" | sed -e s/=/:=/ -e "s/;/:/" > log-td-E0.gappa.autocheck.mpl
echo "printf(\"You should see only zeros in what follows:\\n\");" >> log-td-E0.gappa.autocheck.mpl
sed -f ./TEMPLOG/log-td_1.sed log-td-E0.gappa | $GAPPA 2>&1 | grep simplify >> log-td-E0.gappa.autocheck.mpl
echo "printf(\"If you have seen only zeros up to now, everything's fine \\n\");" >> log-td-E0.gappa.autocheck.mpl

cat log-td-E0-logir0.gappa | grep "#MAPLE" | sed -e s/=/:=/ -e "s/;/:/" > log-td-E0-logir0.gappa.autocheck.mpl
echo "printf(\"You should see only zeros in what follows:\\n\");" >> log-td-E0-logir0.gappa.autocheck.mpl
sed -f ./TEMPLOG/log-td_1.sed log-td-E0-logir0.gappa | $GAPPA 2>&1 | grep simplify >> log-td-E0-logir0.gappa.autocheck.mpl
echo "printf(\"If you have seen only zeros up to now, everything's fine \\n\");" >> log-td-E0-logir0.gappa.autocheck.mpl

cat log-td-accurate.gappa | grep "#MAPLE" | sed -e s/=/:=/ -e "s/;/:/" > log-td-accurate.gappa.autocheck.mpl
echo "printf(\"You should see only zeros in what follows:\\n\");" >> log-td-accurate.gappa.autocheck.mpl
sed -f ./TEMPLOG/log-td-accurate_1.sed log-td-accurate.gappa | $GAPPA 2>&1 | grep simplify >> log-td-accurate.gappa.autocheck.mpl
echo "printf(\"If you have seen only zeros up to now, everything's fine \\n\");" >> log-td-accurate.gappa.autocheck.mpl

cat log-td-accurate-E0.gappa | grep "#MAPLE" | sed -e s/=/:=/ -e "s/;/:/" > log-td-accurate-E0.gappa.autocheck.mpl
echo "printf(\"You should see only zeros in what follows:\\n\");" >> log-td-accurate-E0.gappa.autocheck.mpl
sed -f ./TEMPLOG/log-td-accurate_1.sed log-td-accurate-E0.gappa | $GAPPA 2>&1 | grep simplify >> log-td-accurate-E0.gappa.autocheck.mpl
echo "printf(\"If you have seen only zeros up to now, everything's fine \\n\");" >> log-td-accurate-E0.gappa.autocheck.mpl

cat log-td-accurate-E0-logir0.gappa | grep "#MAPLE" | sed -e s/=/:=/ -e "s/;/:/" > log-td-accurate-E0-logir0.gappa.autocheck.mpl
echo "printf(\"You should see only zeros in what follows:\\n\");" >> log-td-accurate-E0-logir0.gappa.autocheck.mpl
sed -f ./TEMPLOG/log-td-accurate_1.sed log-td-accurate-E0-logir0.gappa | $GAPPA 2>&1 | grep simplify >> log-td-accurate-E0-logir0.gappa.autocheck.mpl
echo "printf(\"If you have seen only zeros up to now, everything's fine \\n\");" >> log-td-accurate-E0-logir0.gappa.autocheck.mpl



echo "read \"log-td.gappa.autocheck.mpl\";" | maple
echo "read \"log-td-E0.gappa.autocheck.mpl\";" | maple
echo "read \"log-td-E0-logir0.gappa.autocheck.mpl\";" | maple
echo "read \"log-td-accurate.gappa.autocheck.mpl\";" | maple
echo "read \"log-td-accurate-E0.gappa.autocheck.mpl\";" | maple
echo "read \"log-td-accurate-E0-logir0.gappa.autocheck.mpl\";" | maple

