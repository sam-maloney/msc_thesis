#!/bin/bash
for limiter in firstOrder; do
#for limiter in firstOrder BarthJespersen; do
#for limiter in firstOrder Venkatakrishnan; do
  for balanced in false; do
    for eps in 0; do
      echo
      echo "Beginning convergence test: ${limiter} ${balanced} ${eps}"
      ./runConvTest.sh ${limiter} ${balanced} ${eps}
      echo
    done
  done
done
