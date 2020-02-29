LOG_FILE=${CASE}.log
PRM_FILE=${CASE}.prm
DAT_FILE=${CASE}.dat

sed -e "s/<np1>/$NP1/g" \
    -e "s/<np2>/$NP2/g" \
    -e "s/<np3>/$NP3/g" \
    -e "s/<ep1>/$EP1/g" \
    -e "s/<ep2>/$EP2/g" \
    -e "s/<ep3>/$EP3/g" \
    -e "s/<po_u>/$PO_U/g" \
    -e "s/<po_p>/$PO_P/g" \
    -e "s/<po_q>/$PO_Q/g" \
    isp_flow__sdc_test.tmpl > ${PRM_FILE}

date > ${LOG_FILE}

srun ../../isp_flow__sdc_test ${CASE} 2>&1 | tee -a ${LOG_FILE}

date >> ${LOG_FILE}

