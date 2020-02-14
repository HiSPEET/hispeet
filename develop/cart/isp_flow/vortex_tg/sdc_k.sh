#!/usr/bin/env bash

# SDC(M,K) with fixed number of subintervals M and varying number of sweeps K

# execute within slurm script or stand-alone using
# `bash -l sdc_m_k.sh`

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

# SDC parameters
N_SUB=${N_SUB:-"3"}
RANGE_N_SWEEP=${RANGE_N_SWEEP:-"0 1 2 3 4 5 6"}

for N_SWEEP in $RANGE_N_SWEEP; do

  CASE=sdc_${N_SUB}-${N_SWEEP}
  LOG_FILE=${CASE}.log
  PRM_FILE=${CASE}.prm
  DAT_FILE=${CASE}.dat

  date > ${LOG_FILE}

  # investigated range of time step sizes
  # dt = 0.1 / sqrt(2^k), k = k_min, ... k_max

  case "$N_SWEEP" in              # trial # study  #
      0)  k_min=3; k_max="25" ;;  #  "8"  #  "25"  #
      1)  k_min=3; k_max="23" ;;  #  "7"  #  "23"  #
      2)  k_min=3; k_max="13" ;;  #  "6"  #  "13"  #
      3)  k_min=3; k_max="11" ;;  #  "5"  #  "11"  #
      4)  k_min=3; k_max="9"  ;;  #  "4"  #  "8"   #
      5)  k_min=3; k_max="7"  ;;  #  "3"  #  "7"   #
      6)  k_min=3; k_max="6"  ;;  #  "3"  #  "6"   #
      *)  k_min=3; k_max="5"      #  "3"  #  "5"   #
  esac

  for((k=k_min; k<=k_max; k++)); do

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
          -e "s/<dt>/$DT/g" \
          -e "s/<n_sub>/$N_SUB/g" \
          -e "s/<n_sweep>/$N_SWEEP/g" \
          isp_flow__sdc_test.tmpl > ${PRM_FILE}

      ${EXEC} ${PROGRAM} ${CASE} 2>&1 | tee -a ${LOG_FILE}

  done

  grep -e "#      t" -m 1 ${LOG_FILE} >  ${DAT_FILE}
  grep -e "#last#$"       ${LOG_FILE} >> ${DAT_FILE}

  date >> ${LOG_FILE}

done

