#!/bin/bash
g++ computeFields.cpp -o computeFields -O3 -std=c++11
echo *x/ ref_12800/ | xargs -n 1 cp computeFields
