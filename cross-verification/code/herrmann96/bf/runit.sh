#!/bin/sh

rm ./*.sac

hprep96 -M modelplate.d -d dfile_BF_R01R03 -HS 1e-6 -HR 99e-6 -TF -BF -ALL -V -NDEC 1
hspec96
hpulse96 -D -IMP -t -l 1 | f96tosac -B

# rename data
rename s/B001/BF_R01_/ ./*.sac
rename s/B002/BF_R03_/ ./*.sac
rm ./*.dat
rm ./*.grn

hprep96 -M modelplate.d -d dfile_BF_R02 -HS 1e-6 -HR 1e-6 -TF -BF -ALL -V -NDEC 1
hspec96
hpulse96 -D -IMP -t -l 1 | f96tosac -B

rename s/B001/BF_R02_/ ./*.sac
