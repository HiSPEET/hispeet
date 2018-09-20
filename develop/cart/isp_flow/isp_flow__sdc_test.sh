#!/usr/bin/env bash

# start with
# > bash -l isp_flow__sdc_test.sh

module load compiler/intel

N_SUB="3"
I_PRC="0"
RANGE_N_CYC="4 5 6"
RANGE_DT=".0009765625 .00048828125 .000244140625 .0001220703125 .00006103515625"

for N_CYC in $RANGE_N_CYC; do
for DT in $RANGE_DT; do

    echo "========================================================================="
    echo "n_cyc = " $N_CYC ", dt =" $DT
    echo
    sed -e "s/<n_sub>/$N_SUB/g" \
        -e "s/<n_cyc>/$N_CYC/g" \
        -e "s/<i_prc>/$I_PRC/g" \
        -e "s/<dt>/$DT/g" \
        isp_flow__sdc_test.tmpl > \
        isp_flow__sdc_test.prm
    #
    mpirun -n 1 isp_flow__sdc_test

done
done
