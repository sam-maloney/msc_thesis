#!/bin/bash
cd $1
#for N in 25 50 100 200 400 800 1600 3200 6400; do
for N in 10 30 90 270 810 2430 7290; do
    delta=`echo "2.0/${N}-2" | bc -l`
    sed -i -e "1 s/.* /$delta /" ${N}_rho
done
cd ..
