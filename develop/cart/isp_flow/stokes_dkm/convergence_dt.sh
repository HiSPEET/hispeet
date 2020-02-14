#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

# number of partitions in directions 1-2
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP=$((${NP1} * ${NP2}))

# use preset execution command, or mpirun, if not set
EXEC=${EXEC:-"mpirun -n ${NP}"}

# path to program
PROGRAM="../isp_flow__sdc_test"

# test case
CASE="convergence_dt"

# elements per partition in directions 1-2 and polynomial orders for u and p
EP=${EP:-"4"}
PO_U="6"
PO_P="5"

# investigated range of time step sizes, dt = 0.1 / 2^k, k = 1, ... 12
RANGE_DT=\
".05 .025 .0125 6.25d-3 3.125d-3 1.5625d-3 7.8125d-4 3.90625d-4 1.953125d-4 9.765625d-5"
#" 4.8828125d-5 2.44140625d-5"

# SDC parameters
N_SUB="1"
RANGE_N_SWEEP="0 1"

for N_SWEEP in $RANGE_N_SWEEP; do

  LOG_FILE=${CASE}-sdc-${N_SUB}-${N_SWEEP}.log
  DAT_FILE=${CASE}-sdc-${N_SUB}-${N_SWEEP}.dat

  date > ${LOG_FILE}

  for DT in $RANGE_DT; do

      echo "========================================================================="
      echo "dt =" $DT
      echo
      sed -e "s/<np1>/$NP1/g" \
          -e "s/<np2>/$NP2/g" \
          -e "s/<ep1>/$EP/g" \
          -e "s/<ep2>/$EP/g" \
          -e "s/<po_u>/$PO_U/g" \
          -e "s/<po_p>/$PO_P/g" \
          -e "s/<dt>/$DT/g" \
          -e "s/<n_sub>/$N_SUB/g" \
          -e "s/<n_sweep>/$N_SWEEP/g" \
          isp_flow__sdc_test.tmpl > \
          isp_flow__sdc_test.prm

      ${EXEC} ${PROGRAM} 2>&1 | tee -a ${LOG_FILE}

  done

  grep -e "#      t" -m 1 ${LOG_FILE} >  ${DAT_FILE}
  grep -e ^" 1.00000E-01" ${LOG_FILE} >> ${DAT_FILE}

  date >> ${LOG_FILE}

done
