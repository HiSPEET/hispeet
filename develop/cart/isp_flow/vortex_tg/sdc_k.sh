#!/usr/bin/env bash

# SDC(M,K) with fixed number of subintervals M and varying number of sweeps K

# execute within slurm script or stand-alone using
# `bash -l sdc_k.sh`

# number of partitions in directions 1-2
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP=$((${NP1} * ${NP2}))

# use preset execution command, or mpirun, if not set
MPIRUN=${MPIRUN:-"mpirun"}
EXEC=${EXEC:-"${MPIRUN} -n ${NP}"}

# path to program
PROGRAM="../isp_flow__sdc_test"

# elements per partition in directions 1-2
EP=${EP:-"2"}

# polynomial orders for u and p and for integrating the nonlinear terms
PO_U="16"
PO_P="15"
PO_Q="24"

# final time
T_END=${T_END:-"0.03125"}

# SDC parameters
N_SUB=${N_SUB:-"3"}
N_SWEEP_MAX=$((3*${N_SUB}))

# max time step size and max number of time step subdivisions per series
DT_MAX=${DT_MAX:-"0.015625"}
ST_MAX=${ST_MAX:-"4"}

for((N_SWEEP=0; N_SWEEP<=N_SWEEP_MAX; N_SWEEP++)); do

  CASE=sdc_${N_SUB}-${N_SWEEP}
  LOG_FILE=${CASE}.log
  PRM_FILE=${CASE}.prm
  DAT_FILE=${CASE}.dat

  date > ${LOG_FILE}

  # investigated range of time step sizes
  # dt = DT_MAX / sqrt(2^s), s = 0, ... ST

  case "$N_SWEEP" in            
      0)  ST="10" ;;
      1)  ST="10" ;;
      2)  ST="10" ;;
      3)  ST="10" ;;
      4)  ST="10" ;;
      5)  ST="10" ;;
      6)  ST="9"  ;;
      *)  ST="8"    
  esac

  ST=$((ST > ST_MAX ? ST_MAX : ST))

  for((s=0; s<=ST; s++)); do

      DT=$(bc -l <<< "${DT_MAX}/sqrt(2^${s})")

      echo "========================================================================="
      echo "M =" ${N_SUB} ", K =" ${N_SWEEP} ", dt =" $DT, ", s =" $s
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

