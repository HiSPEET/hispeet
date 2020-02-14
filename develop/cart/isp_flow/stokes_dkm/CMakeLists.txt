file( COPY isp_flow__sdc_test.tmpl
           stokes_dkm.prm
           convergence_dt.sh
           convergence_dt.py
           stability_dt.sh
           stability_inf-sup.sh
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (SLURM_FOUND)

  if (OpenMP)
      set( RUN_CONVERGENCE_DT    convergence_dt__slurm_omp.sh    )
      set( RUN_STABILITY_DT      stability_dt__slurm_omp.sh      )
      set( RUN_STABILITY_INF-SUP stability_inf-sup__slurm_omp.sh )
  elseif (OpenACC)
      set( RUN_CONVERGENCE_DT    convergence_dt__slurm_acc.sh    )
      set( RUN_STABILITY_DT      stability_dt__slurm_acc.sh      )
      set( RUN_STABILITY_INF-SUP stability_inf-sup__slurm_acc.sh )
  else ()
      set( RUN_CONVERGENCE_DT    convergence_dt__slurm.sh        )
      set( RUN_STABILITY_DT      stability_dt__slurm.sh          )
      set( RUN_STABILITY_INF-SUP stability_inf-sup__slurm.sh     )
  endif ()

  configure_file( ${RUN_CONVERGENCE_DT}    convergence_dt__slurm.sh    )
  configure_file( ${RUN_STABILITY_DT}      stability_dt__slurm.sh      )
  configure_file( ${RUN_STABILITY_INF-SUP} stability_inf-sup__slurm.sh )

endif (SLURM_FOUND)
