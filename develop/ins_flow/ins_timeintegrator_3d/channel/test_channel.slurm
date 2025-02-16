#!/bin/bash

# ==============================================================================
# SLURM script for channel example using MPI
# ------------------------------------------------------------------------------

##SBATCH --account=                          # project
#SBATCH --nodes=2
#SBATCH --ntasks=128                         # max number of MPI-processes
#SBATCH --time=00-16:00:00                   # dd-hh:mm:ss
#SBATCH --job-name=channel                   # sensible job name

# ------------------------------------------------------------------------------

if [ X"$SLURM_STEP_ID" = "X" -a X"$SLURM_PROCID" = "X"0 ]
then
  echo ""
  echo "SLURM_JOB_ID        = $SLURM_JOB_ID"
  echo "SLURM_JOB_PARTITION = $SLURM_JOB_PARTITION"
  echo "SLURM_JOB_NODELIST  = $SLURM_JOB_NODELIST"
  echo ""
fi

CASE=test_1g
#CASE=test_mg

# ------------------------------------------------------------------------------
# Execution

date

module load hispeet/gnu12

if ! compgen -G "*.h5" > /dev/null; then
    echo "starting ==========================================================================="
    cp ${CASE}_start.prm ${CASE}.prm
    srun --ntasks $SLURM_NTASKS --cpu_bind=cores --distribution=block:block \
         ../ins_timeintegrator_3d_test ${CASE}
    mkdir -p save
    cp *.h5 save/
fi

for ((k = 1; k <= 10; k++)); do
    echo "restarting, k =" $k "=================================================================="
    cp ${CASE}_restart.prm ${CASE}.prm
    srun --ntasks $SLURM_NTASKS --cpu_bind=cores --distribution=block:block \
         ../ins_timeintegrator_3d_test ${CASE}
    cp *h5 save/
done

date


