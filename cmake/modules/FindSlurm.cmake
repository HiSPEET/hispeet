find_program(Slurm_SBATCH_COMMAND sbatch DOC "Path to the Slurm sbatch executable")
find_program(Slurm_SRUN_COMMAND   srun   DOC "Path to the Slurm srun executable")

find_package_handle_standard_args( Slurm DEFAULT_MSG 
                                   Slurm_SBATCH_COMMAND 
                                   Slurm_SRUN_COMMAND   )

mark_as_advanced(Slurm_SRUN_COMMAND Slurm_SBATCH_COMMAND)
