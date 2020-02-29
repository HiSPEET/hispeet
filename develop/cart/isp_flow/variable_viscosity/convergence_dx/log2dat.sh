#!/bin/bash

for PO in "06" "11" "16"
do
    echo "P =" ${PO}
    DAT_FILE="p${PO}.dat"
    rm -f ${DAT_FILE}

    for LOG_FILE in p${PO}*.log
    do
        echo $LOG_FILE  "-->" $DAT_FILE
        if [ ! -f ${DAT_FILE} ]
        then
            grep -e "#      t" -m 1 ${LOG_FILE} >  ${DAT_FILE}
        fi
        grep -e "#last#$" ${LOG_FILE} >> ${DAT_FILE}
    done
done

