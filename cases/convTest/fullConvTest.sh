#!/bin/bash
for limiter in firstOrder BarthJespersen; do
#for limiter in firstOrder Venkatakrishnan; do
  for balanced in true; do
    for eps in 0 0.0001 0.1 2; do
      echo
      echo "Beginning convergence test: ${limiter} ${balanced} ${eps}"
      ./runConvTest.sh ${limiter} ${balanced} ${eps}
      echo
    done
  done
done
