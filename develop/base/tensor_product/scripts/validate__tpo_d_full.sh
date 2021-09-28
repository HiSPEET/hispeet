##!/bin/bash

# ==============================================================================
# Run script for evaluating a 1:1 tensor-product operator
#

# ------------------------------------------------------------------------------
# Settings

TEST_CASES="grad_3d_d"

# Set NE_TEST to perform tests of NP_RANGE
NE_TEST=8
NP_RANGE="3 4 5"

# Set NP_TEST to perform tests of NE_RANGE
NP_TEST=12
NE_RANGE="1 2 3 4"

# Set number of tests x operand size for adapting number of test runs
NT_NO=10000000

# min/max number of test runs
NT_MIN=1
NT_MAX=1

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

  # test over range of elements for each operator size..........................

  if [ -n "$NP_TEST" ]
  then

    echo ""
    echo "test over range of elements for each operator size"

    FILE=${PROGRAM}"_np"$NP".dat"

    for NP in $NP_RANGE; do

      (( PO = NP - 1  ))
      echo "  np =" $NP
      for NE in $NE_RANGE; do
        echo "     nr/nz =" $NE

      # set number of test runs
      NT=$(( NT_NO / (NP*NP*NP*NE) ))
      NT=$(( NT > NT_MAX ? NT_MAX : NT ))
      NT=$(( NT < NT_MIN ? NT_MIN : NT ))

      sed -e "s/<po>/$PO/g" \
          -e "s/<nr>/$NE/g" \
          -e "s/<np>/$NE/g" \
          -e "s/<nz>/$NE/g" \
          -e "s/<nt>/$NT/g" \
            ${PROGRAM}.tmpl > ${PROGRAM}.prm

      $EXEC $PROGRAM | grep ^[[:blank:]]*[1-9] >> ${FILE}

      done
    done
  fi


done

echo ""
echo "done!"
echo ""
