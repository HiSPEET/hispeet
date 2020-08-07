#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

# test case
CASE="convergence_dt"

# path to program
PROGRAM="../isp_flow__sdc_test"

# number of partitions in directions 1-2
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP=$((${NP1} * ${NP2}))

# use preset execution command, or mpirun, if not set
MPIRUN=${MPIRUN:-"mpirun"}
EXEC=${EXEC:-"${MPIRUN} -n ${NP}"}

# elements per partition in directions 1-2 and polynomial orders for u and p
EP=${EP:-"4"}
PO_U="10"
PO_P="9"

# time integration
T_END=${T_END:-"0.25"}
TIME_METHOD=${TIME_METHOD:-"5"}
N_SUB_MAX=${N_SUB_MAX:-"2"}

# max time step size and max number of time step subdivisions per series
DT_MAX=${DT_MAX:-"0.025"}
ST_MAX=${ST_MAX:-"6"}

for((N_SUB=0; N_SUB<=N_SUB_MAX; N_SUB++)); do

  N_SWEEP=$((3 * ${N_SUB}))

  LOG_FILE=${CASE}-sdc-${N_SUB}.log
  DAT_FILE=${CASE}-sdc-${N_SUB}.dat

  date > ${LOG_FILE}

  for((s=0; s<=ST_MAX; s++)); do
  
      DT=$(bc -l <<< "${DT_MAX}/sqrt(2^${s})")
  
      echo "========================================================================="
      echo "dt =" $DT ", s = "${s}"/"${ST_MAX}
      echo
      sed -e "s/<np1>/$NP1/g" \
          -e "s/<np2>/$NP2/g" \
          -e "s/<ep1>/$EP/g" \
          -e "s/<ep2>/$EP/g" \
          -e "s/<po_u>/$PO_U/g" \
          -e "s/<po_p>/$PO_P/g" \
          -e "s/<dt>/$DT/g" \
          -e "s/<time_method>/$TIME_METHOD/g" \
          -e "s/<n_sub>/$N_SUB/g" \
          -e "s/<n_sweep>/$N_SWEEP/g" \
          isp_flow__sdc_test.tmpl > \
          isp_flow__sdc_test.prm

      ${EXEC} ${PROGRAM} 2>&1 | tee -a ${LOG_FILE}

  done

  grep -e "#      t" -m 1 ${LOG_FILE} >  ${DAT_FILE}
  grep -e "#last#"        ${LOG_FILE} >> ${DAT_FILE}

  date >> ${LOG_FILE}

done
