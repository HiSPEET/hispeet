#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

PROG="../ins_timeintegrator_3d_test"
EXEC=${EXEC:-"mpirun"}
SERIES=${SERIES:-"stokes_dkm"}

# range for time steps DT = 1 / 10^k
KMIN=${KMIN:-"1"}
KMAX=${KMAX:-"5"}

# number of stages of substeps
NS=${NS:-"1"}

# reference velocity
V_REF="3.46307"

# reference bulk viscosity
MU_REF=${MU_REF:-"1"}

# method for computing MU
#
#   1)  MU = MU_REF
#   2)  MU = MU_REF * V_REF * DELTA
#   3)  MU = MU_REF * V_REF * DELTA * max(1, 1/CS)
#
# where 
#   DELTA = DX/2 * DELTA_S 
#   CS    = V_REF * DELTA / (DT * NS)
# with NS stages or substeps
#
MU_METHOD=${MU_METHOD:-"0"}

# polynomial orders
PO_U=${PO_U:-"4"}
PO_P=$((${PO_U} - 1))

# convective eigenvalue λ for standard element of order P = PO_U
case $PO_U in
   2)  LAMBDA="1.00000000000";;
   3)  LAMBDA="1.66481452034";;
   4)  LAMBDA="2.47159803598";;
   5)  LAMBDA="3.40883748057";;
   6)  LAMBDA="4.47073758438";;
   7)  LAMBDA="5.65599485263";;
   8)  LAMBDA="6.96732548733";;
   9)  LAMBDA="8.41089833351";;
  10)  LAMBDA="9.99507344979";;
  11)  LAMBDA="11.7283555837";;
  12)  LAMBDA="13.6174001000";;
  13)  LAMBDA="15.6661755899";;
  14)  LAMBDA="17.8764016028";;
  15)  LAMBDA="20.2484698510";;
  16)  LAMBDA="22.7821848530";;
  17)  LAMBDA="25.4771710732";;
  18)  LAMBDA="28.3330470912";;
  19)  LAMBDA="31.3494832099";;
  20)  LAMBDA="34.5262112697";;
  21)  LAMBDA="37.8630182223";;
  22)  LAMBDA="41.3597359477";;
  23)  LAMBDA="45.0162316204";;
  24)  LAMBDA="48.8323997653";;
  25)  LAMBDA="52.8081560476";;
  26)  LAMBDA="56.9434325203";;
  27)  LAMBDA="61.2381740125";;
  28)  LAMBDA="65.6923353835";;
  29)  LAMBDA="70.3058794266";;
  30)  LAMBDA="75.0787752594";;
  31)  LAMBDA="80.0109970785";;
  32)  LAMBDA="85.1025231895";;
   *)  echo "PO_U out of range";;
esac

# test grids
if [[ $HPC ]]
then
    # for HPC system
    NP=( 1 1 2 4 8 )  # number of partitions in directions 1-2 
    EP=( 1 2 2 2 2 )  # number of elements per partition in directions 1-2 
else
    # for personal computers or laptops, using 4 cores
    NP=( 1 1 2 2 2 )  # number of partitions in directions 1-2 
    EP=( 1 2 2 4 8 )  # number of elements per partition in directions 1-2 
fi
NG=${#NP[@]}  # number of grids used

# times step sizes DT = 1 / 10^k, K_MIN ≤ k ≤ K_MAX
for k in $(seq $KMIN $KMAX) ; do

    DT=$(bc -l <<< "1 / 10^$k")

    # grids 0 .. NG-1 
    for i in $(seq 0 $(($NG - 1))) ; do
	
        CASE=${SERIES}"_po"${PO_U}"_i"${i}"_k"${k}
        LOGFILE=${CASE}".log"
        
        date > ${LOGFILE}
        
        NE=$(bc -l <<< "${NP[$i]} * ${EP[$i]}")
        DX=$(bc -l <<< "2/${NE}")
        DELTA=$(bc -l <<< "${DX} / (2 * ${LAMBDA})")
        CFL=$(bc -l <<< "${V_REF} * ${DT} / ${DELTA}")
        CS=$(bc -l <<< "${CFL}/${NS}")

        case $MU_METHOD in
          1)  MU=${MU_REF};;
          2)  MU=$(bc -l <<< "${MU_REF} * ${V_REF} * ${DELTA}");;
          3)  MU=$(bc -l <<< "${MU_REF} * ${V_REF} * ${DELTA}")
              MU=$(bc -l <<< "if (${CS} < 1) ${MU}/${CS} else ${MU}");;
        esac

        NPROC=$(bc -l <<< "${NP[$i]} * ${NP[$i]}")         

        echo                   >> ${LOGFILE}
        echo "NE    =" $NE     >> ${LOGFILE}
        echo "DX    =" $DX     >> ${LOGFILE}
        echo "DELTA =" $DELTA  >> ${LOGFILE}
        echo "DT    =" $DT     >> ${LOGFILE}
        echo "CFL   =" $CFL    >> ${LOGFILE}
        echo "CS    =" $CS     >> ${LOGFILE}
        echo "MU    =" $MU     >> ${LOGFILE}
        echo "NPROC =" $NPROC  >> ${LOGFILE}

        sed -e "s/<dx>/$DX/g" \
            -e "s/<mu>/$MU/g" \
            -e "s/<np>/${NP[$i]}/g" \
            -e "s/<ep>/${EP[$i]}/g" \
            -e "s/<po_u>/$PO_U/g" \
            -e "s/<po_p>/$PO_P/g" \
            -e "s/<dt>/$DT/g" \
            study.tmpl > ${CASE}.prm

        ${EXEC} -n ${NPROC} ${PROG} ${CASE} 2>&1 | tee -a ${LOGFILE}

        echo >> ${LOGFILE}
        date >> ${LOGFILE}

    done
done
