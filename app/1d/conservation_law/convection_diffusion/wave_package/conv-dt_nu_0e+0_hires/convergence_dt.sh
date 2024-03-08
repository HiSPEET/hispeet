#!/usr/bin/env bash

# test case
CASE="convergence_dt"

# path to program
PROGRAM="../../../conservation_law_1l"

# diffusivity
NU="0"

# final time
T_END="10.0"

# CFL range 2^K_MIN ... 2^K_MAX
K_MIN="-6" 
#K_MAX="8"
K_MAX="0"

# time integration method
# 1  IMEX Euler
# 2  SI1 (ISD1)
# 3  SI2 (ISD2)
# 4  IMEX Runge-Kutta of order 3, 5 stages (ARS443)
# 5  TVD Runge-Kutta of order 3
TIME_METHOD="2"

# SDC method
# 0  none
# 1  IMEX Euler
# 2  SI1
SDC_METHOD="2"

# SDC point set, num collocation points, num correction sweeps, num stages
SDC_POINT_SET="RR"
SDC_N_COL="8"
SDC_N_SWEEP="14"
SDC_N_STAGE="1"

# time step size for CFL = 2^0
DT_0="3.8583162370057488532e-4" 

rm -f ${CASE}.err

date > ${CASE}.log

for((k = K_MIN; k <= K_MAX; k++)); do

    CFL=$(bc -l <<< "1 / 2^${k}")
    DT=$(bc -l <<< "${DT_0} * ${CFL}")

    echo "========================================================================="
    echo "k =" $k", CFL =" $CFL", dt =" $DT
    echo
    sed -e "s/<NU>/$NU/g" \
        -e "s/<T_END>/$T_END/g" \
        -e "s/<DT>/$DT/g" \
        -e "s/<TIME_METHOD>/$TIME_METHOD/g" \
        -e "s/<SDC_METHOD>/$SDC_METHOD/g" \
        -e "s/<SDC_POINT_SET>/$SDC_POINT_SET/g" \
        -e "s/<SDC_N_COL>/$SDC_N_COL/g" \
        -e "s/<SDC_N_SWEEP>/$SDC_N_SWEEP/g" \
        -e "s/<SDC_N_STAGE>/$SDC_N_STAGE/g" \
        ${CASE}.tmpl > \
        ${CASE}.prm

    ${PROGRAM} ${CASE} 2>&1 | tee -a ${CASE}.log
    
    if [ -f "${CASE}.run" ]; then
        if [ -f "${CASE}.err" ]; then
            tail -n 1 ${CASE}.run >> ${CASE}.err
            rm ${CASE}.run
        else
            mv ${CASE}.run ${CASE}.err
        fi
    fi

done

date >> ${CASE}.log

