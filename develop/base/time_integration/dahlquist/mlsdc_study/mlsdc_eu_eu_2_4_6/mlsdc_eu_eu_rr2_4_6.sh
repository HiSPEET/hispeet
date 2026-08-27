#!/usr/bin/env bash

# test case
CASE="mlsdc_eu_eu_rr2_4_6"

# path to program
PROGRAM="../../dahlquist_mlsdc"

# parameters

P_TIME1="2"
P_TIME2="4"
P_TIME3="6"
MG_START=$1
N_COARSE=$2
C_MAX="14"        # CYCLE_MAX = 2*N_COL2 - 2
PRE_METHOD="1"
PRE_ORDER="1"
PRE_STAGES="1"
SDC_METHOD="1"
SDC_POINTS="'RR'"
SDC_STAGES="1"

# test range for amplification
C_MAX_AMP="20"
D_MIN_AMP="-5"
D_MAX_AMP="10"

# test range for accuracy
C_MAX_ACC=".6"
D_MIN_ACC="-.50"
D_MAX_ACC=".75"

# test dimensions
NC="2000"
ND="2000"

date > ${CASE}.log

# amplification
for((C = 0; C <= C_MAX; C++)); do

    echo "========================================================================="
    echo "amplification: C =" $C
    echo
    sed -e "s/<PRE_METHOD>/$PRE_METHOD/g" \
        -e "s/<PRE_ORDER>/$PRE_ORDER/g" \
        -e "s/<PRE_STAGES>/$PRE_STAGES/g" \
        -e "s/<SDC_METHOD>/$SDC_METHOD/g" \
        -e "s/<SDC_POINTS>/$SDC_POINTS/g" \
        -e "s/<SDC_STAGES>/$SDC_STAGES/g" \
        -e "s/<MG_START>/$MG_START/g" \
        -e "s/<N_COARSE>/$N_COARSE/g" \
        -e "s/<P_TIME1>/$P_TIME1/g" \
	    -e "s/<P_TIME2>/$P_TIME2/g" \
        -e "s/<P_TIME3>/$P_TIME3/g" \
        -e "s/<CYC>/$C/g" \
        -e "s/<C_MAX>/$C_MAX_AMP/g" \
        -e "s/<D_MIN>/$D_MIN_AMP/g" \
        -e "s/<D_MAX>/$D_MAX_AMP/g" \
        -e "s/<NC>/$NC/g" \
        -e "s/<ND>/$ND/g" \
        ${CASE}.tmpl > \
	${CASE}.prm

    ${PROGRAM} ${CASE} >&1 | tee -a ${CASE}.log
    
    mv lambda_im.dat lambda_im_amp.dat
    mv lambda_re.dat lambda_re_amp.dat
    mv lambda_im.dat lambda_im_acc.dat
    mv lambda_re.dat lambda_re_acc.dat
    mv amplification_level_3.dat amplification_level_3  ${MG_START}_cyc${C}.dat
    mv error_level_3.dat error_level_3_acc_mg${MG_START}_cyc${C}.dat

done

date >> ${CASE}.log
