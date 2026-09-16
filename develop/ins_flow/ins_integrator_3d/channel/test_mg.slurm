#!/bin/bash

# ==============================================================================
# SLURM script for channel example using MPI
# ------------------------------------------------------------------------------

#SBATCH --account=p_hispeet
#SBATCH --nodes=1
#SBATCH --ntasks=64                          # max number of MPI-processes
#SBATCH --time=00-01:00:00                   # dd-hh:mm:ss
#SBATCH --job-name=channel_mg                # sensible job name

# ------------------------------------------------------------------------------

if [ X"$SLURM_STEP_ID" = "X" -a X"$SLURM_PROCID" = "X"0 ]
then
  echo ""
  echo "SLURM_JOB_ID        = $SLURM_JOB_ID"
  echo "SLURM_JOB_PARTITION = $SLURM_JOB_PARTITION"
  echo "SLURM_JOB_NODELIST  = $SLURM_JOB_NODELIST"
  echo ""
fi

PROG=${PROG:-"../ins_integrator_3d_test"}
CASE=${CASE:-"test_mg"}

RESTART_TAG_IN=${RESTART_TAG_IN:-init}
RESTART_TAG_OUT=${RESTART_TAG_OUT:-save}

# number of restart loops
N_RESTART=${N_RESTART:-"1"}

# averaging rate (0 disabled)
AVG_RATE=${AVG_RATE:-"0"}

# max number of iterations in diffusion solver
I_MAX_V=${I_MAX_V:-"50"}

# step size
DT=${DT:-"2e-5"}

# max number of steps
NT_MAX=${NT_MAX:-"1000"}

echo "PROG            =" $PROG
echo "CASE            =" $CASE
echo "N_RESTART       =" $N_RESTART
echo "RESTART_TAG_IN  =" $RESTART_TAG_IN
echo "RESTART_TAG_OUT =" $RESTART_TAG_OUT

# ------------------------------------------------------------------------------
# Execution

date

#module load hispeet/barnard/gnu12
#module load hispeet/barnard/gnu14

if ! compgen -G "*${RESTART_TAG_IN}*.h5" > /dev/null; then
    echo "starting ==========================================================================="
    # create prm-file
    sed -e "s/<restart_tag_in>//g" \
        -e "s/<restart_tag_out>/$RESTART_TAG_OUT/g" \
        -e "s/<avg_rate>/$AVG_RATE/g" \
        -e "s/<i_max_v>/$I_MAX_V/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<nt_max>/$NT_MAX/g" \
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
        -e "s/<avg_rate>/$AVG_RATE/g" \
        -e "s/<i_max_v>/$I_MAX_V/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<nt_max>/$NT_MAX/g" \
        ${CASE}.tmpl > ${CASE}.prm
    # run program
    srun --ntasks $SLURM_NTASKS --cpu_bind=cores --distribution=block:block \
         ${PROG} ${CASE}
done

date

