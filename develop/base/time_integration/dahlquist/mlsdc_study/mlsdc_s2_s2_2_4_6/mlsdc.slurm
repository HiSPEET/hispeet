#!/bin/bash

#SBATCH --job-name=mlsdc_study
#SBATCH --exclusive
#SBATCH --output=out_slurm_%j.log
#SBATCH --error=error_slurm_%j.log
#SBATCH --mail-type=all
#SBATCH --time=10:00:00                 #h:m:s
#SBATCH --account=p_hispeet

./mlsdc_study.sh
