file( COPY isp_flow__sdc_test.tmpl
           vortex_hw.prm
           convergence_dt.sh
           convergence_dt.py
           convergence_dx.sh
           convergence_dx.py
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (Slurm_FOUND)

  if (OpenMP)
      set( RUN_CONVERGENCE_DT convergence_dt__omp.slurm )
      set( RUN_CONVERGENCE_DX convergence_dx__omp.slurm )
  elseif (OpenACC)
      set( RUN_CONVERGENCE_DT convergence_dt__acc.slurm )
      set( RUN_CONVERGENCE_DX convergence_dx__acc.slurm )
  else ()
      set( RUN_CONVERGENCE_DT convergence_dt.slurm )
      set( RUN_CONVERGENCE_DX convergence_dx.slurm )
  endif ()

  configure_file( ${RUN_CONVERGENCE_DT} convergence_dt.slurm )
  configure_file( ${RUN_CONVERGENCE_DX} convergence_dx.slurm )

endif (Slurm_FOUND)
