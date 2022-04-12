#!/bin/bash

# ==============================================================================
# Run slurm scripts with different number of threads

# ------------------------------------------------------------------------------
# Settings

THREAD_RANGE="1 2 3 4 5 6 7 8 9 10 11 12"
DEPENDENCY=""
JOB_FILE="validate__tpo_d.slurm" # "validate__tpo_r.slurm"

for THREADS in $THREAD_RANGE ; do
    JOB_CMD="sbatch --cpus-per-task=$THREADS"
    if [ -n "$DEPENDENCY" ] ; then
        JOB_CMD="$JOB_CMD --dependency afterany:$DEPENDENCY"
    fi
    JOB_CMD="$JOB_CMD $JOB_FILE"
    echo -n "Running command: $JOB_CMD  "
    OUT=`$JOB_CMD`
    echo "Result: $OUT"
    DEPENDENCY=`echo $OUT | awk '{print $4}'`
done
