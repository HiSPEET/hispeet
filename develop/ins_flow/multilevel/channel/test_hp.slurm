#!/bin/bash

# ==============================================================================
# SLURM script for channel example using MPI
# ------------------------------------------------------------------------------

#SBATCH --account=p_hispeet
#SBATCH --nodes=1
#SBATCH --ntasks=64                          # max number of MPI-processes
#SBATCH --time=00-05:00:00                   # dd-hh:mm:ss
#SBATCH --job-name=channel_ml                # sensible job name

# ------------------------------------------------------------------------------

if [ X"$SLURM_STEP_ID" = "X" -a X"$SLURM_PROCID" = "X"0 ]
then
  echo ""
  echo "SLURM_JOB_ID        = $SLURM_JOB_ID"
  echo "SLURM_JOB_PARTITION = $SLURM_JOB_PARTITION"
  echo "SLURM_JOB_NODELIST  = $SLURM_JOB_NODELIST"
  echo ""
fi

PROG=${PROG:-"../ml_ins_solver_3d"}
CASE=${CASE:-"test_hp"}


RESTART_TAG_IN=${RESTART_TAG_IN:-init}
RESTART_TAG_OUT=${RESTART_TAG_OUT:-save}
RESTART_ADJUST=${RESTART_ADJUST:-"F"}

# number of restart loops
N_RESTART=${N_RESTART:-"1"}

# averaging rate (0 disabled)
AVG_RATE=${AVG_RATE:-"0"}

# initial top level
L_TOP=${L_TOP:-"5"}

# switch to mixed order
MIXED=${MIXED:-"F"}

# Grad-Siv stabilization
DIV_STAB=${DIV_STAB:-"1"}
MU_0=${MU_0:-"1"}

# SGS model
SGS_MODEL=${SGS_MODEL:-"0"}
SGS_DYNAMIC=${SGS_DYNAMIC:-"F"}

# step size, may be increased in sucessive runs
DT=${DT:-"1e-5"}

# max number of steps
NT_MAX=${NT_MAX:-"1000"}


echo "PROG            =" $PROG
echo "CASE            =" $CASE
echo "RESTART_TAG_IN  =" $RESTART_TAG_IN
echo "RESTART_TAG_OUT =" $RESTART_TAG_OUT
echo "RESTART_ADJUST  =" $RESTART_ADJUST
echo "N_RESTART       =" $N_RESTART
echo "L_TOP           =" $L_TOP
echo "MIXED           =" $MIXED
echo "DIV_STAB        =" $DIV_STAB
echo "MU_0            =" $MU_0
echo "SGS_MODEL       =" $SGS_MODEL
echo "SGS_DYNAMIC     =" $SGS_DYNAMIC
echo "DT              =" $DT
echo "NT_MAX          =" $NT_MAX

# ------------------------------------------------------------------------------
# Execution

date

if ! compgen -G "*${RESTART_TAG_IN}*.h5" > /dev/null; then
    echo "starting ==========================================================================="
    # create prm-file
    sed -e "s/<restart_tag_in>//g" \
        -e "s/<restart_tag_out>/$RESTART_TAG_OUT/g" \
        -e "s/<restart_adjust>/$RESTART_ADJUST/g" \
        -e "s/<avg_rate>/$AVG_RATE/g" \
        -e "s/<l_top>/$L_TOP/g" \
        -e "s/<mixed>/$MIXED/g" \
        -e "s/<div_stab>/$DIV_STAB/g" \
        -e "s/<mu_0>/$MU_0/g" \
        -e "s/<sgs_model>/$SGS_MODEL/g" \
        -e "s/<sgs_dynamic>/$SGS_DYNAMIC/g" \
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
    for outfile in *${RESTART_TAG_OUT}*.*; do
        infile=${outfile/$RESTART_TAG_OUT/$RESTART_TAG_IN}
        cp $outfile $infile
	    mv $outfile save/
    done
    # creating prm-file
    sed -e "s/<restart_tag_in>/$RESTART_TAG_IN/g" \
        -e "s/<restart_tag_out>/$RESTART_TAG_OUT/g" \
        -e "s/<restart_adjust>/$RESTART_ADJUST/g" \
        -e "s/<avg_rate>/$AVG_RATE/g" \
        -e "s/<l_top>/$L_TOP/g" \
        -e "s/<mixed>/$MIXED/g" \
        -e "s/<div_stab>/$DIV_STAB/g" \
        -e "s/<mu_0>/$MU_0/g" \
        -e "s/<sgs_model>/$SGS_MODEL/g" \
        -e "s/<sgs_dynamic>/$SGS_DYNAMIC/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<nt_max>/$NT_MAX/g" \
        ${CASE}.tmpl > ${CASE}.prm
    # run program
    srun --ntasks $SLURM_NTASKS --cpu_bind=cores --distribution=block:block \
         ${PROG} ${CASE}
done

date

