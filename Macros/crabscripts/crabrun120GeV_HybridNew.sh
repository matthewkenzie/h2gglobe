#!/bin/sh

mkdir outputToy
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
./combine 120GeVmodel_3sigma.root -M HybridNew -T $n  --generateNuis=0 --generateExt=1 --fitNuisances=1 --testStat Atlas -D data_mass -m 120 -s -1 -t 1 --generateBinnedWorkaround -S 1 -H ProfileLikelihood --hintStatOnly
rm CMS-HGG.root
rm 120GeVmodel_3sigma.root
mv *.root outputToy/
rm *.root
tar cvfz outputToy.tgz outputToy/
