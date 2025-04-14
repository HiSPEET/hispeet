#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

CASE="stokes_dkm"
PROG="../ins_timeintegrator_3d_test"
EXEC=${EXEC:-"mpirun"}
NPROC_MAX=${NPROC_MAX:-"4"}

# polynomial orders and time step
PO_U=${PO_U:-"4"}
PO_P=$((${PO_U} - 1))
DT=${DT:-"1e-4"}

# Peclet number for stabilization
PE_MU=${PE_MU:-"1"}

# bulk viscosity μ = vδ with v = v_ref and δ(Δx=1,P=po_u)
case $PO_U in
   2)  MU_REF="1.7315349578857422";;
   3)  MU_REF="1.0400767993869648";;
   4)  MU_REF="0.70057304330241035";;
   5)  MU_REF="0.50795468183945190";;
   6)  MU_REF="0.38730409137342146";;
   7)  MU_REF="0.30614153707719081";;
   8)  MU_REF="0.24852218559818247";;
   9)  MU_REF="0.20586801661690601";;
  10)  MU_REF="0.17323884277430450";;
  11)  MU_REF="0.14763663546306542";;
  12)  MU_REF="0.12715606100790033";;
  13)  MU_REF="0.11052697245396699";;
  14)  MU_REF="0.096861493513128083";;
  15)  MU_REF="0.085514360869113305";;
  16)  MU_REF="0.076003902569398019";;
   *)  echo "PO_U out of range";;
esac

echo "========================================================================="
echo
date
echo

NP=1
EP=2
NE=$((${NP}*${EP}))
L3=$(bc -l <<< "2/${NE}")
MU=$(bc -l <<< "${MU_REF}/${PE_MU} * 2/${NE}")

echo "ne =" $NE
echo "dx =" $L3
echo "mu =" $MU

sed -e "s/<l3>/$L3/g" \
    -e "s/<mu>/$MU/g" \
    -e "s/<np>/${NP}/g" \
    -e "s/<ep>/${EP}/g" \
    -e "s/<po_u>/$PO_U/g" \
    -e "s/<po_p>/$PO_P/g" \
    -e "s/<dt>/$DT/g" \
    test_convergence.tmpl > \
    test_convergence.prm

${EXEC} -n 1 ${PROG} test_convergence 2>&1 | tee ${CASE}"_p"${PO_U}"_n"${NE}.log

echo "========================================================================="
echo
date
echo

NP=2
EP=2
NE=$((${NP}*${EP}))
L3=$(bc -l <<< "2/${NE}")
MU=$(bc -l <<< "${MU_REF}/${PE_MU} * 2/${NE}")

echo "ne =" $NE
echo "dx =" $L3
echo "mu =" $MU

sed -e "s/<l3>/$L3/g" \
    -e "s/<mu>/$MU/g" \
    -e "s/<np>/${NP}/g" \
    -e "s/<ep>/${EP}/g" \
    -e "s/<po_u>/$PO_U/g" \
    -e "s/<po_p>/$PO_P/g" \
    -e "s/<dt>/$DT/g" \
    test_convergence.tmpl > \
    test_convergence.prm

${EXEC} -n 4 ${PROG} test_convergence 2>&1 | tee ${CASE}"_p"${PO_U}"_n"${NE}.log

echo "========================================================================="
echo
date
echo

if [ "$NPROC_MAX" -ge "16" ]
then
  NP=4
  EP=2
else
  NP=2
  EP=4
fi
NPROC=$((${NP}*${NP}))

NE=$((${NP}*${EP}))
L3=$(bc -l <<< "2/${NE}")
MU=$(bc -l <<< "${MU_REF}/${PE_MU} * 2/${NE}")

echo "ne =" $NE
echo "dx =" $L3
echo "mu =" $MU

sed -e "s/<l3>/$L3/g" \
    -e "s/<mu>/$MU/g" \
    -e "s/<np>/${NP}/g" \
    -e "s/<ep>/${EP}/g" \
    -e "s/<po_u>/$PO_U/g" \
    -e "s/<po_p>/$PO_P/g" \
    -e "s/<dt>/$DT/g" \
    test_convergence.tmpl > \
    test_convergence.prm

${EXEC} -n ${NPROC} ${PROG} test_convergence 2>&1 | tee ${CASE}"_p"${PO_U}"_n"${NE}.log

echo "========================================================================="
echo
date
echo

if [ "$NPROC_MAX" -ge "64" ]
then
  NP=8
  EP=2
else
  NP=2
  EP=8
fi
NPROC=$((${NP}*${NP}))

NE=$((${NP}*${EP}))
L3=$(bc -l <<< "2/${NE}")
MU=$(bc -l <<< "${MU_REF}/${PE_MU} * 2/${NE}")

echo "ne =" $NE
echo "dx =" $L3
echo "mu =" $MU

sed -e "s/<l3>/$L3/g" \
    -e "s/<mu>/$MU/g" \
    -e "s/<np>/${NP}/g" \
    -e "s/<ep>/${EP}/g" \
    -e "s/<po_u>/$PO_U/g" \
    -e "s/<po_p>/$PO_P/g" \
    -e "s/<dt>/$DT/g" \
    test_convergence.tmpl > \
    test_convergence.prm

${EXEC} -n ${NPROC} ${PROG} test_convergence 2>&1 | tee ${CASE}"_p"${PO_U}"_n"${NE}.log
