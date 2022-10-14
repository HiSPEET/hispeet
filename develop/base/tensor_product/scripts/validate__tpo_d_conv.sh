#!/bin/bash

# ==============================================================================
# Run script for evaluating the incompressible Navier-Stokes convection operator

# ------------------------------------------------------------------------------
# Settings

TEST_CASES="ins_convection_3d_d"

# Set polynomial order of mesh
PO_G=3

# Set NE_TEST to perform tests of NP_RANGE
NE_TEST=8
NV_RANGE="8 12"

# Set NP_TEST to perform tests of NE_RANGE
NV_TEST=
NE_RANGE="3 5"

# Set number of tests x operand size for adapting number of test runs
NT_NO=10000000

# min/max number of test runs
NT_MIN=10
NT_MAX=10

# ------------------------------------------------------------------------------
# Execution

for CASE in $TEST_CASES; do

  PROGRAM="./validate__tpo__"${CASE}

  if [ ! -f $PROGRAM ]; then
    echo $PROGRAM "does not exist -- skipping test"
    continue
  fi

  echo ""
  echo "----------------------------------------------------------------"
  echo "Evaluation of tensor-product operator "${CASE}

  # test over range of elements ................................................

  if [ -n "$NV_TEST" ]
  then

    (( NV = NV_TEST ))
    (( PO_V = NV - 1  ))

    echo ""
    echo "test over range of elements"
      sed -e "s/<po_v>/$PO_V/g" \
    (( NQ = ((3 * NV) + 1) / 2 ))
    (( PO_Q = NQ - 1))

    FILE=${PROGRAM}"_nv"$NV".dat"

    for NE in $NE_RANGE; do

      echo "  ne =" $NE

      # set number of test runs
      NT=$(( NT_NO / (NV*NV*NV*NE) ))
      NT=$(( NT > NT_MAX ? NT_MAX : NT ))
      NT=$(( NT < NT_MIN ? NT_MIN : NT ))

      sed -e "s/<pog>/$PO_G/g" \
          -e "s/<pov>/$PO_V/g" \
          -e "s/<poq>/$PO_Q/g" \
          -e "s/<nr>/$NE/g" \
          -e "s/<np>/$NE/g" \
          -e "s/<nz>/$NE/g" \
          -e "s/<nt>/$NT/g" \
            ${PROGRAM}.tmpl > ${PROGRAM}.prm

      $EXEC $PROGRAM | grep ^[[:blank:]]*[1-9] >> ${FILE}

    done
  fi

  # test over range of operator sizes ..........................................

  if [ -n "$NE_TEST" ]
  then

    echo ""
    echo "test over range of operator sizes"

    (( NE = NE_TEST ))

    FILE=${PROGRAM}"_ne"$NE".dat"

    for NV in $NV_RANGE; do

        (( PO_V = NV - 1 ))
        (( NQ = ((3 * NV) + 1) / 2 ))
        (( PO_Q = NQ - 1))
        echo "  po_v =" $PO_V



        # set number of test runs
        NT=$(( NT_NO / (NV*NV*NV*NE) ))
        NT=$(( NT > NT_MAX ? NT_MAX : NT ))
        NT=$(( NT < NT_MIN ? NT_MIN : NT ))

      sed -e "s/<pog>/$PO_G/g" \
          -e "s/<pov>/$PO_V/g" \
          -e "s/<poq>/$PO_Q/g" \
          -e "s/<nr>/$NE/g" \
          -e "s/<np>/$NE/g" \
          -e "s/<nz>/$NE/g" \
          -e "s/<nt>/$NT/g" \
            ${PROGRAM}.tmpl > ${PROGRAM}.prm

      $EXEC $PROGRAM | grep ^[[:blank:]]*[1-9] >> ${FILE}

    done
  fi

done

echo ""
echo "done!"
echo ""
