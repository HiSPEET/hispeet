#!/usr/bin/env bash

# path to program
PROGRAM="./convection_diffusion"

# test case
CASE="convergence_dt"

# time integration method
METHOD="6"

# range of sweeps with SDC
RANGE_N_SWEEP="0 1 2"

#  investigated range of time steps: dt = ct / 2^k, k = k_min, .. k_max
ct="0.01"
k_min="2"
k_max="5"

for N_SWEEP in $RANGE_N_SWEEP; do

    SUBCASE=${CASE}"_"${N_SWEEP}

    date > ${SUBCASE}.log
    
    for ((k=k_min; k<=k_max; k++)); do
    
        DT=$(bc -l <<< "${ct} / 2^${k}")
    
        sed -e "s/<method>/$METHOD/g" \
            -e "s/<dt>/$DT/g" \
            -e "s/<n_sweep>/$N_SWEEP/g" \
            convergence_dt.tmpl > \
            convergence_dt.prm
    
        ${PROGRAM} convergence_dt 2>&1 | tee -a ${SUBCASE}.log
    
    done
    date >> ${SUBCASE}.log
 
    grep -e "#      t" -m 1 ${SUBCASE}.log >  ${SUBCASE}.dat
    grep -e "#last#$"       ${SUBCASE}.log >> ${SUBCASE}.dat

done
