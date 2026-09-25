#!/bin/bash

# ==============================================================================
# SLURM script for test_scalar_product
# ------------------------------------------------------------------------------

#SBATCH --account=p_hispeet          # project
#SBATCH --nodes=8                    # number of nodes
#SBATCH --ntasks-per-node=104        # max number of MPI processes
##SBATCH --cpus-per-task=1           # number of cores per process with OpenMP
##SBATCH --cpu-freq=highm1           # CPU frequency in kHz
##SBATCH --mem=60000                 # memory per node in MB

#SBATCH --time=00-05:00:00
##SBATCH --exclusive
#SBATCH --job-name=test_scalar_product      
#SBATCH --output=test_scalar_product_%j.out

# ------------------------------------------------------------------------------

if [ X"$SLURM_STEP_ID" = "X" -a X"$SLURM_PROCID" = "X"0 ]
then
  echo ""
  echo "SLURM_JOB_ID        = $SLURM_JOB_ID"
  echo "SLURM_JOB_PARTITION = $SLURM_JOB_PARTITION"
  echo "SLURM_JOB_NODELIST  = $SLURM_JOB_NODELIST"
  echo ""
fi

# ------------------------------------------------------------------------------
# Scaling tests
#
# The number of MPI processes starts at 1 and is doubled until the maximum
# (nodes x tasks per node) is reached. If the maximum is not a power of two,
# a final run with the maximum is appended. Each run uses the smallest number
# of nodes that can host its processes.
#
#   Part 1 - strong scaling: fixed total vector length N_STRONG
#   Part 2 - weak scaling:   fixed vector length per process N_WEAK,
#                            i.e. total length n_proc * N_WEAK

N_STRONG=120000000                     # total vector length (strong scaling)
N_WEAK=200000                          # vector length per process (weak scaling)
N_TEST=100                             # repetitions per run (mean time)

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
N_NODES=${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-1}}
TASKS_PER_NODE=${SLURM_NTASKS_PER_NODE:-$(( ${SLURM_NTASKS:-1} / N_NODES ))}
MAX_PROCS=$(( N_NODES * TASKS_PER_NODE ))

# process counts: 1, 2, 4, ... (, MAX_PROCS)
PROC_LIST=()
for (( np = 1; np <= MAX_PROCS; np *= 2 )); do PROC_LIST+=($np); done
[ ${PROC_LIST[-1]} -ne $MAX_PROCS ] && PROC_LIST+=($MAX_PROCS)

echo "Nodes = $N_NODES, tasks/node = $TASKS_PER_NODE, max procs = $MAX_PROCS," \
     "threads/proc = $OMP_NUM_THREADS"
echo "Process counts: ${PROC_LIST[*]}"

# ------------------------------------------------------------------------------
# run_series MODE RESULTS
#   MODE = strong : n = N_STRONG,         speedup = t(1)/t(p), effic = speedup/p
#   MODE = weak   : n = n_proc * N_WEAK,  effic = t(1)/t(p)

run_series() {
  local mode=$1 results=$2 t1="" np nn n OUT t_ref t_hyb

  if [ $mode = strong ]; then
    printf '# %8s %6s %8s %12s %13s %13s %9s %9s\n' \
           n_proc nodes threads n t_ref t_hybrid speedup effic > $results
  else
    printf '# %8s %6s %8s %12s %13s %13s %9s\n' \
           n_proc nodes threads n t_ref t_hybrid effic > $results
  fi

  for np in "${PROC_LIST[@]}"; do

    nn=$(( (np + TASKS_PER_NODE - 1) / TASKS_PER_NODE ))  # nodes needed
    if [ $mode = strong ]; then n=$N_STRONG; else n=$(( np * N_WEAK )); fi

    echo ""
    echo "========== $mode: n_proc = $np on $nn node(s), n = $n =========="

    OUT=$(srun --nodes=$nn --ntasks=$np --cpus-per-task=$OMP_NUM_THREADS \
               ./test_scalar_product $n $N_TEST)
    echo "$OUT"

    t_ref=$(echo "$OUT" | awk '/reference:/ {print $NF}')
    t_hyb=$(echo "$OUT" | awk '/hybrid:/    {print $NF}')
    [ -z "$t_hyb" ] && { echo "run with $np processes failed"; continue; }
    [ -z "$t1" ] && t1=$t_hyb

    awk -v mode=$mode -v np=$np -v nn=$nn -v nt=$OMP_NUM_THREADS -v n=$n \
        -v tr=$t_ref -v th=$t_hyb -v t1=$t1 'BEGIN {
          if (mode == "strong") {
            s = t1 / th
            printf "  %8d %6d %8d %12d %13.5e %13.5e %9.2f %9.3f\n", \
                   np, nn, nt, n, tr, th, s, s/np
          } else {
            printf "  %8d %6d %8d %12d %13.5e %13.5e %9.3f\n", \
                   np, nn, nt, n, tr, th, t1/th
          }
        }' >> $results
  done
}

# ------------------------------------------------------------------------------
# Part 1: strong scaling

STRONG_RESULTS=strong_scaling_${SLURM_JOB_ID:-local}.dat
echo ""
echo "#################### Part 1: strong scaling, n = $N_STRONG ####################"
run_series strong $STRONG_RESULTS

# ------------------------------------------------------------------------------
# Part 2: weak scaling

WEAK_RESULTS=weak_scaling_${SLURM_JOB_ID:-local}.dat
echo ""
echo "#################### Part 2: weak scaling, n / n_proc = $N_WEAK ####################"
run_series weak $WEAK_RESULTS

# ------------------------------------------------------------------------------
# Summary

echo ""
echo "==================== Strong scaling summary ===================="
cat $STRONG_RESULTS
echo ""
echo "==================== Weak scaling summary ======================"
cat $WEAK_RESULTS
