#!/bin/sh

rm ./*.sac

hprep96 -M modelplate.d -d dfile_DC_R01R03 -HS 50e-6 -HR 99e-6 -TF -BF -ALL -V -NDEC 1
hspec96
hpulse96 -V -STEP -t -l 1 | f96tosac -B

# rename data
rename s/B001/DC_R01_/ ./*.sac
rename s/B002/DC_R03_/ ./*.sac
rm ./*.dat
rm ./*.grn

hprep96 -M modelplate.d -d dfile_DC_R02 -HS 50e-6 -HR 1e-6 -TF -BF -ALL -V -NDEC 1
hspec96
hpulse96 -V -STEP -t -l 1 | f96tosac -B

rename s/B001/DC_R02_/ ./*.sac
