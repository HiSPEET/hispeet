file( COPY isp_flow__sdc_test.tmpl
           isp_flow__sdc_test.prm
           vortex_tg.prm
           sdc_k.sh
           sdc_m.sh
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (SLURM_FOUND)

  if (OpenMP)
      set( RUN_SDC_K sdc_k__slurm_omp.sh )
      set( RUN_SDC_M sdc_m__slurm_omp.sh )
  elseif (OpenACC)
      set( RUN_SDC_K sdc_k__slurm_acc.sh )
      set( RUN_SDC_M sdc_m__slurm_acc.sh )
  else ()
      set( RUN_SDC_K sdc_k__slurm.sh )
      set( RUN_SDC_M sdc_m__slurm.sh )
  endif ()

  configure_file( ${RUN_SDC_K} sdc_k__slurm.sh )
  configure_file( ${RUN_SDC_M} sdc_m__slurm.sh )

endif (SLURM_FOUND)
