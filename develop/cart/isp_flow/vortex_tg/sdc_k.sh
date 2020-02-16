#!/usr/bin/env bash

# SDC(M,K) with fixed number of subintervals M and varying number of sweeps K

# execute within slurm script or stand-alone using
# `bash -l sdc_k.sh`

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
N_SUB=${N_SUB:-"3"}
N_SWEEP_MAX=$((3*${N_SUB}))

for((N_SWEEP=0; N_SWEEP<=N_SWEEP_MAX; N_SWEEP++)); do

  CASE=sdc_${N_SUB}-${N_SWEEP}
  LOG_FILE=${CASE}.log
  PRM_FILE=${CASE}.prm
  DAT_FILE=${CASE}.dat

  date > ${LOG_FILE}

  # investigated range of time step sizes
  # dt = 0.1 / sqrt(2^k), k = k_min, ... k_max

  case "$N_SWEEP" in            
      0)  k_min=2; k_max="12" ;;
      1)  k_min=2; k_max="12" ;;
      2)  k_min=2; k_max="12" ;;
      3)  k_min=2; k_max="10" ;;
      4)  k_min=2; k_max="9"  ;;
      5)  k_min=2; k_max="8"  ;;
      6)  k_min=2; k_max="7"  ;;
      *)  k_min=2; k_max="6"    
  esac

  for((k=k_min; k<=k_max; k++)); do

      DT=$(bc -l <<< "0.1/sqrt(2^${k})")

      echo "========================================================================="
      echo "M =" ${N_SUB} ", K =" ${N_SWEEP} ", dt =" $DT
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
          isp_flow__sdc_test.tmpl > ${PRM_FILE}

      ${EXEC} ${PROGRAM} ${CASE} 2>&1 | tee -a ${LOG_FILE}

  done

  grep -e "#      t" -m 1 ${LOG_FILE} >  ${DAT_FILE}
  grep -e "#last#$"       ${LOG_FILE} >> ${DAT_FILE}

  date >> ${LOG_FILE}

done

