file( COPY README.md
           convergence_dt.py
           convergence_dt.sh
           convergence_dx.py
           convergence_dx.sh
           ins_timeintegrator_3d_test.prm
           ins_timeintegrator_3d_test.tmpl
           vortex_tg.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)

  if (OpenMP)
      configure_file(convergence_dt-hyb.slurm convergence_dt.slurm)
      configure_file(convergence_dx-hyb.slurm convergence_dx.slurm)
  else ()
      configure_file(convergence_dt-mpi.slurm convergence_dt.slurm)
      configure_file(convergence_dx-mpi.slurm convergence_dx.slurm)
  endif ()

endif (Slurm_FOUND)
