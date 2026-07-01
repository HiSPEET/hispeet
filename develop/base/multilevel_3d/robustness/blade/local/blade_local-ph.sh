#!/usr/bin/env bash

# ==============================================================================
# Run script for blade multigrid test

# ------------------------------------------------------------------------------
# Settings

# executable
PROGRAM="../../../ml_elliptic_test_static"

# test case
SERIES="blade_local-ph"

# mpi
EXEC=${EXEC:-"mpirun"}
NPROC=${NPROC:-"8"}

# parameters
NPR=${NPR:-"8"}
NPG=${NPG:-"8"}
IMAX=${IMAX:-"10"}
NS1=${NS1:-"6"}
NS2=${NS2:-"2"}
NSC=${NSC:-"4"}

# case name
CASE=${SERIES}"-ns_"${NS1}"_"${NS2}"_"${NSC}

# ------------------------------------------------------------------------------
# Execution

printf "\n"
printf "=%.0s" {1..100}
printf "\n%s%s\n\n" "Executing case " ${CASE}

sed -e "s/<npr>/$NPR/g" \
    -e "s/<npg>/$NPG/g" \
    -e "s/<imax>/$IMAX/g" \
    -e "s/<ns1>/$NS1/g" \
    -e "s/<ns2>/$NS2/g" \
    -e "s/<nsc>/$NSC/g" \
    ${SERIES}.tmpl > ${CASE}.prm

$EXEC -n $NPROC $PROGRAM ${CASE} 2>&1 | tee ${CASE}-np_${NPROC}.log

echo ""
echo "done!"
echo ""
