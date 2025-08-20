#!/bin/bash

echo '"piS","piT","k","m","replicate"' > analysis/pi.csv
echo '"count","k","m","replicate"' > analysis/sfs.csv

cat simulation/data/pi/*.csv >> analysis/pi.csv
cat simulation/data/sfs/*.csv >> analysis/sfs.csv