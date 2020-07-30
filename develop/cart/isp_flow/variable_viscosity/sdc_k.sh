#!/usr/bin/env bash

# SDC(M,K) with fixed number of subintervals M and varying number of sweeps K

# execute within slurm script or stand-alone using
# `bash -l sdc_k.sh`

# path to program
PROGRAM="../isp_flow__sdc_test"

# number of partitions in directions 1-3
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP3=${NP3:-"1"}
NP=$((${NP1} * ${NP2} * ${NP3}))

# use preset execution command, or mpirun, if not set
EXEC=${EXEC:-"mpirun -n ${NP}"}

# elements per partition in directions 1-3
EP1=${EP1:-"2"}
EP2=${EP2:-"2"}
EP3=${EP3:-"2"}

# polynomial orders for u and p and for integrating the nonlinear terms
PO_U="16"
PO_P="15"
PO_Q="24"

# time integration
T_END=${T_END:-"0.25"}
TIME_METHOD=${TIME_METHOD:-"5"}

# SDC parameters
N_SUB=${N_SUB:-"3"}
RANGE_N_SWEEP=${RANGE_N_SWEEP:-"0 1 2 3 4 5 6"}

# investigated range of time step sizes, dt = dt_0 / 2^k, k = k_min, ... k_max
dt_0="0.1"
k_min="2"
k_max="6"

for N_SWEEP in $RANGE_N_SWEEP; do

  CASE=sdc_${N_SUB}-${N_SWEEP}
  LOG_FILE=${CASE}.log
  PRM_FILE=${CASE}.prm
  DAT_FILE=${CASE}.dat

  date > ${LOG_FILE}

  for ((k=k_min; k<=k_max; k++)); do

      DT=$(bc -l <<< "${dt_0}/sqrt(2^${k})")

      echo "========================================================================="
      echo "M =" ${N_SUB} ", K =", ${N_SWEEP} ", dt =" $DTT  "(",$k,"/",$k_max,")"
      echo
      echo
      sed -e "s/<np1>/$NP1/g" \
          -e "s/<np2>/$NP2/g" \
          -e "s/<np3>/$NP3/g" \
          -e "s/<ep1>/$EP1/g" \
          -e "s/<ep2>/$EP2/g" \
          -e "s/<ep3>/$EP3/g" \
          -e "s/<po_u>/$PO_U/g" \
          -e "s/<po_p>/$PO_P/g" \
          -e "s/<po_q>/$PO_Q/g" \
          -e "s/<t_end>/$T_END/g" \
          -e "s/<dt>/$DT/g" \
          -e "s/<time_method>/$TIME_METHOD/g" \
          -e "s/<n_sub>/$N_SUB/g" \
          -e "s/<n_sweep>/$N_SWEEP/g" \
          isp_flow__sdc_test.tmpl > ${PRM_FILE}

        ${EXEC} ${PROGRAM} ${CASE} 2>&1 | tee -a ${LOG_FILE}

  done

  grep -e "#      t" -m 1 ${LOG_FILE} >  ${DAT_FILE}
  grep -e "#last#$"       ${LOG_FILE} >> ${DAT_FILE}

  date >> ${LOG_FILE}

done

