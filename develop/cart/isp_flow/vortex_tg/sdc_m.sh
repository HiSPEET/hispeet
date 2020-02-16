#!/usr/bin/env bash

# SDC(M,K=3*M) with varying number of subintervals M

# execute within slurm script or stand-alone using
# `bash -l sdc_m.sh`

# number of partitions in directions 1-2
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP=$((${NP1} * ${NP2}))

# use preset execution command, or mpirun, if not set
EXEC=${EXEC:-"mpirun -n ${NP}"}

# path to program
PROGRAM="../isp_flow__sdc_test"

# elements per partition in directions 1-2
EP=${EP:-"2"}

# polynomial orders for u and p and for integrating the nonlinear terms
PO_U="16"
PO_P="15"
PO_Q="24"

# Final time
T_END=${T_END:-"0.1"}

# SDC parameters
N_SUB_MAX=${N_SUB_MAX:-"2"}

for((N_SUB=0; N_SUB<=N_SUB_MAX; N_SUB++)); do

  N_SWEEP=$((3*${N_SUB}))

  LOG_FILE=sdc_${N_SUB}-${N_SWEEP}.log
  DAT_FILE=sdc_${N_SUB}-${N_SWEEP}.dat

  date > ${LOG_FILE}

  # investigated range of time step sizes
  # dt = 0.1 / sqrt(2^k), k = 2, ... k_max

  case "$N_SUB" int
      1)  k_max="12" ;;
      2)  k_max="10" ;;
      3)  k_max="8"  ;;
      4)  k_max="6"  ;;
      *)  k_max="4"
  esac

  for((k=2; k<=k_max; k++)); do

      DT=$(bc -l <<< "0.1/sqrt(2^${k})")

      echo "========================================================================="
      echo "M =" ${N_SUB} ", K =", ${N_SWEEP} ", dt =" $DT
      echo
      sed -e "s/<np1>/$NP1/g" \
          -e "s/<np2>/$NP2/g" \
          -e "s/<ep1>/$EP/g" \
          -e "s/<ep2>/$EP/g" \
          -e "s/<po_u>/$PO_U/g" \
          -e "s/<po_p>/$PO_P/g" \
          -e "s/<po_q>/$PO_Q/g" \
          -e "s/<t_end>/$T_END/g" \
          -e "s/<dt>/$DT/g" \
          -e "s/<n_sub>/$N_SUB/g" \
          -e "s/<n_sweep>/$N_SWEEP/g" \
          isp_flow__sdc_test.tmpl > \
          isp_flow__sdc_test.prm

      ${EXEC} ${PROGRAM} 2>&1 | tee -a ${LOG_FILE}

  done

  grep -e "#      t" -m 1 ${LOG_FILE} >  ${DAT_FILE}
  grep -e "#last#$"       ${LOG_FILE} >> ${DAT_FILE}

  date >> ${LOG_FILE}

done
