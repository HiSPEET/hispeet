#!/bin/bash

# ==============================================================================
# Run script for evaluating a 1:1 tensor-product operator
#

# ------------------------------------------------------------------------------
# Settings

TEST_CASES="grad_3d_r"

# Set NE_TEST to perform tests of NP_RANGE
NE_TEST=1024
NP_RANGE="8 12 16"

# Set NP_TEST to perform tests of NE_RANGE
NP_TEST=16
NE_RANGE="8 16 24 32 48 64 128 256 512 1024 2048"

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

  if [ -n "$NP_TEST" ]
  then

    echo ""
    echo "test over range of elements"

    (( NP = NP_TEST ))
    (( PO = NP - 1  ))

    FILE=${PROGRAM}"_np"$NP".dat"

    for NE in $NE_RANGE; do

      echo "  ne =" $NE

      # set number of test runs
      NT=$(( NT_NO / (NP*NP*NP*NE) ))
      NT=$(( NT > NT_MAX ? NT_MAX : NT ))
      NT=$(( NT < NT_MIN ? NT_MIN : NT ))

      sed -e "s/<po>/$PO/g" \
          -e "s/<np>/$NP/g" \
          -e "s/<ne>/$NE/g" \
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

    for NP in $NP_RANGE; do

        echo "  po =" $PO

        (( PO = NP - 1 ))

        # set number of test runs
        NT=$(( NT_NO / (NP*NP*NP*NE) ))
        NT=$(( NT > NT_MAX ? NT_MAX : NT ))
        NT=$(( NT < NT_MIN ? NT_MIN : NT ))

      sed -e "s/<po>/$PO/g" \
          -e "s/<np>/$NP/g" \
          -e "s/<ne>/$NE/g" \
          -e "s/<nt>/$NT/g" \
            ${PROGRAM}.tmpl > ${PROGRAM}.prm

      $EXEC $PROGRAM | grep ^[[:blank:]]*[1-9] >> ${FILE}

    done
  fi

done

echo ""
echo "done!"
echo ""
