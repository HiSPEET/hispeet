file( COPY isp_flow__sdc_test.tmpl
           stokes_dkm.prm
           convergence_dt.sh
           convergence_dt.py
           stability_dt.sh
           stability_dt.py
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (SLURM_FOUND)

  if (OpenMP)
      set( RUN_CONVERGENCE_DT  convergence_dt__omp.slurm )
      set( RUN_STABILITY_DT    stability_dt__omp.slurm   )
  elseif (OpenACC)
      set( RUN_CONVERGENCE_DT  convergence_dt__acc.slurm )
      set( RUN_STABILITY_DT    stability_dt__acc.slurm   )
  else ()
      set( RUN_CONVERGENCE_DT  convergence_dt.slurm      )
      set( RUN_STABILITY_DT    stability_dt.slurm        )
  endif ()

  configure_file( ${RUN_CONVERGENCE_DT}    convergence_dt__slurm.sh    )
  configure_file( ${RUN_STABILITY_DT}      stability_dt__slurm.sh      )

endif (SLURM_FOUND)
