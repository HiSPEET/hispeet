#!/usr/bin/env bash

# ==============================================================================
# Run script for blade multigrid test

# ------------------------------------------------------------------------------
# Settings

# executable
PROGRAM="../../../ml_elliptic_test_static"

# test case
SERIES="blade_global-ph"


# mpi
EXEC=${EXEC:-"mpirun"}
echo $NPROC
NPROC=${NPROC:-"8"}
echo $NPROC

# parameters
NPR=${NPR:-"8"}
NPG=${NPG:-"8"}
START=${START:-"3"}
ICRS=${ICRS:-"10000"}
IMAX=${IMAX:-"10"}
NS1=${NS1:-"6"}
NS2=${NS2:-"2"}
NSC=${NSC:-"4"}
RCRS=${RCRS:-"1e-10"}
RRED=${RRED:-"1e-10"}
RMAX=${RMAX:-"1e-12"}

# case name
CASE=${SERIES}"-ns_"${NS1}"_"${NS2}"_"${NSC}

# ------------------------------------------------------------------------------
# Execution

printf "\n"
printf "=%.0s" {1..100}
printf "\n%s%s\n\n" "Executing case " ${CASE}

sed -e "s/<npr>/$NPR/g" \
    -e "s/<npg>/$NPG/g" \
    -e "s/<start>/$START/g" \
    -e "s/<icrs>/$ICRS/g" \
    -e "s/<imax>/$IMAX/g" \
    -e "s/<ns1>/$NS1/g" \
    -e "s/<ns2>/$NS2/g" \
    -e "s/<nsc>/$NSC/g" \
    -e "s/<rcrs>/$RCRS/g" \
    -e "s/<rred>/$RRED/g" \
    -e "s/<rmax>/$RMAX/g" \
    ${SERIES}.tmpl > ${CASE}.prm

$EXEC -n $NPROC $PROGRAM ${CASE} 2>&1 | tee ${CASE}-np_${NPROC}.log

echo ""
echo "done!"
echo ""
