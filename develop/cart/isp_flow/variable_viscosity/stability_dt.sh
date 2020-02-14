#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l stability_dt.sh`

# number of partitions in directions 1-3
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP3=${NP3:-"1"}
NP=$((${NP1} * ${NP2} * ${NP3}))

# use preset execution command, or mpirun, if not set
EXEC=${EXEC:-"mpirun -n ${NP}"}

# path to program
PROGRAM="../isp_flow__sdc_test"

# elements per partition in directions 1-2
EP1=${EP1:-"1"}
EP2=${EP2:-"1"}
EP3=${EP3:-"1"}

# polynomial orders for u and p and for integrating the nonlinear terms
PO_U=${PO_U:-"8"}
PO_P=${PO_P:-"7"}
PO_Q=${PO_Q:-"12"}

# SDC parameters
N_SUB=${N_SUB:-"1"}
N_SWEEP=${N_SWEEP:-"0"}

# time stepping
DT0=${DT0:-"0"}
DT1=${DT1:-"1"}
NBI=${NBI:-"15"}

# test case
CASE="stability_dt__P"${PO_U}"__M"${N_SUB}"__K"${N_SWEEP}

echo "========================================================================="
echo ${CASE}
echo


date > ${CASE}.log


for((i=1; i<=NBI; i++)); do

    DT=$(bc -l <<< "(${DT0} + ${DT1}) / 2")

    echo "dt =" $DT
    sed -e "s/<np1>/$NP1/g" \
        -e "s/<np2>/$NP2/g" \
        -e "s/<np3>/$NP3/g" \
        -e "s/<ep1>/$EP1/g" \
        -e "s/<ep2>/$EP2/g" \
        -e "s/<ep3>/$EP3/g" \
        -e "s/<po_u>/$PO_U/g" \
        -e "s/<po_p>/$PO_P/g" \
        -e "s/<po_q>/$PO_Q/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<n_sub>/$N_SUB/g" \
        -e "s/<n_sweep>/$N_SWEEP/g" \
        isp_flow__sdc_test.tmpl > \
        ${CASE}.prm

    ${EXEC} ${PROGRAM} ${CASE} 2>&1 | tee ${CASE}.run

    # give sytem time for writing to disk
    sleep 60
    
    if grep -q -e failed -e NaN -e Infinity ${CASE}.run; then
        DT1=$DT
    else
        DT0=$DT
    fi

    cat ${CASE}.run >> ${CASE}.log

done

date >> ${CASE}.log
