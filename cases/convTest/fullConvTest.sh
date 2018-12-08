#!/bin/bash
#for limiter in firstOrder BarthJespersen; do
for limiter in firstOrder Venkatakrishnan; do
  for balanced in false true; do
    for eps in 0.1 0.0001 2; do
      echo
      echo "Beginning convergence test: ${limiter} ${balanced} ${eps}"
      ./runConvTest.sh ${limiter} ${balanced} ${eps}
      echo
    done
  done
done
