file( COPY isp_flow__sdc_test.prm
           isp_flow__sdc_test.tmpl
           variable_viscosity.prm
           convergence_dt.sh
           convergence_dx.sh
           stability_dt.sh
           sdc_k.sh
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (SLURM_FOUND)

  if (OpenMP)
      set( RUN_CONVERGENCE_DT convergence_dt__omp.slurm )
      set( RUN_CONVERGENCE_DX convergence_dx__omp.slurm )
      set( RUN_STABILITY_DT   stability_dt__omp.slurm   )
  elseif (OpenACC)
      set( RUN_CONVERGENCE_DT convergence_dt__acc.slurm )
      set( RUN_CONVERGENCE_DX convergence_dx__acc.slurm )
      set( RUN_STABILITY_DT   stability_dt__acc.slurm   )
  else ()
      set( RUN_CONVERGENCE_DT convergence_dt.slurm )
      set( RUN_CONVERGENCE_DX convergence_dx.slurm )
      set( RUN_STABILITY_DT   stability_dt.slurm   )
      set( RUN_SDC_K          sdc_k.slurm          )
  endif ()

  configure_file( ${RUN_CONVERGENCE_DT} convergence_dt.slurm )
  configure_file( ${RUN_CONVERGENCE_DX} convergence_dx.slurm )
  configure_file( ${RUN_STABILITY_DT}   stability_dt.slurm   )
  configure_file( ${RUN_SDC_K}          sdc_k.slurm          )

endif (SLURM_FOUND)


add_subdirectory(convergence_dx)