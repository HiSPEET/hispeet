#!/bin/bash

#SBATCH -J MLSDC_Parameter
#SBATCH --exclusive
#SBATCH --output=out_slurm_%j.log
#SBATCH --error=error_slurm_%j.log
#SBATCH --mail-type=all
#SBATCH --time=10:00:00                 #h:m:s
#SBATCH --account=p_hispeet

module load hispeet/gnu-foss-hdf5-1.14.4

./MLSDC_Parameterstudy.sh
