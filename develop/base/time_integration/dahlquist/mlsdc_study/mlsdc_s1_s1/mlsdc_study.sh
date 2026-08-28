#!/usr/bin/env bash

# Define the parameter ranges
MG_START_VALUES=(1)
N_COARSE_VALUES=(2) # (2 3)

# Iterate over the parameter ranges
for MG_START in ${MG_START_VALUES[@]}; do
  for N_COARSE in ${N_COARSE_VALUES[@]}; do
        # Prepare the script file with the current parameters
        SCRIPT_FILE="mlsdc_s1_s1_rr2_4_6.sh"
        if [[ ! -f $SCRIPT_FILE ]]; then
          echo "Error: $SCRIPT_FILE not found."
          exit 1
        fi

        # Trigger the script with the current parameters and redirect output to log file
        LOG_FILE="mg_start_${MG_START}_n_coarse_${N_COARSE_VALUES}.log"
        ./$SCRIPT_FILE $MG_START $N_COARSE > $LOG_FILE
  done
done
