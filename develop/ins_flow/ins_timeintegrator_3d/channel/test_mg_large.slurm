#!/bin/bash

# ==============================================================================
# SLURM script for channel example using MPI
# ------------------------------------------------------------------------------

#SBATCH --account=p_hispeet
#SBATCH --nodes=5
#SBATCH --ntasks=256                         # max number of MPI-processes
#SBATCH --time=00-01:00:00                   # dd-hh:mm:ss
#SBATCH --job-name=channel_mg_large          # sensible job name

# ------------------------------------------------------------------------------

if [ X"$SLURM_STEP_ID" = "X" -a X"$SLURM_PROCID" = "X"0 ]
then
  echo ""
  echo "SLURM_JOB_ID        = $SLURM_JOB_ID"
  echo "SLURM_JOB_PARTITION = $SLURM_JOB_PARTITION"
  echo "SLURM_JOB_NODELIST  = $SLURM_JOB_NODELIST"
  echo ""
fi

PROG=${PROG:-"../ins_timeintegrator_3d_test"}
CASE=${CASE:-"test_mg_large"}

N_RESTART=${N_RESTART:-1}
RESTART_TAG_IN=${RESTART_TAG_IN:-init}
RESTART_TAG_OUT=${RESTART_TAG_OUT:-save}

echo "PROG            =" $PROG
echo "CASE            =" $CASE
echo "N_RESTART       =" $N_RESTART
echo "RESTART_TAG_IN  =" $RESTART_TAG_IN
echo "RESTART_TAG_OUT =" $RESTART_TAG_OUT

# ------------------------------------------------------------------------------
# Execution

date

#module load hispeet/barnard/gnu12

if ! compgen -G "*${RESTART_TAG_IN}*.h5" > /dev/null; then
    echo "starting ==========================================================================="
    # create prm-file
    sed -e "s/<restart_tag_in>//g" \
        -e "s/<restart_tag_out>/$RESTART_TAG_OUT/g" \
        ${CASE}.tmpl > ${CASE}.prm
    # run program
    srun --ntasks $SLURM_NTASKS --cpu_bind=cores --distribution=block:block \
         ${PROG} ${CASE}
fi

for ((k = 1; k <= ${N_RESTART}; k++)); do
    echo "restarting, k =" $k "=================================================================="
    # renaming + saving restart files
    mkdir -p save
    for outfile in *${RESTART_TAG_OUT}*.h5; do
        infile=${outfile/$RESTART_TAG_OUT/$RESTART_TAG_IN}
        cp $outfile $infile
	    mv $outfile save/
    done
    # creating prm-file
    sed -e "s/<restart_tag_in>/$RESTART_TAG_IN/g" \
        -e "s/<restart_tag_out>/$RESTART_TAG_OUT/g" \
        ${CASE}.tmpl > ${CASE}.prm
    # run program
    srun --ntasks $SLURM_NTASKS --cpu_bind=cores --distribution=block:block \
         ${PROG} ${CASE}
done

date

