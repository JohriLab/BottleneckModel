#! /bin/bash
sed -i '' '1d' analysis/sfs_3.csv
sed -i '' '1d' analysis/pi_3.csv

cat analysis/sfs_3.csv >> analysis/sfs.csv
cat analysis/pi_3.csv >> analysis/pi.csv