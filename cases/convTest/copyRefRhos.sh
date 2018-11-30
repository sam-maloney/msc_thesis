#!/bin/bash
mv ref_rho0 ref_rho
echo [BfV]*e0/ | xargs -n 1 cp ref_rho
mv ref_rho ref_rho0
mv ref_rho0.1 ref_rho
echo [BfV]*e0.1/ | xargs -n 1 cp ref_rho
mv ref_rho ref_rho0.1
mv ref_rho1 ref_rho
echo [BfV]*e1/ | xargs -n 1 cp ref_rho
mv ref_rho ref_rho1
