file( COPY isp_flow__sdc_test.tmpl
           isp_flow__sdc_test.prm
           transition_tg.prm
           sdc_k.sh
           sdc_k.py
           sdc_m.sh
           sdc_m.py
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (SLURM_FOUND)

  if (OpenMP)
      set( RUN_SDC_K sdc_k__omp.slurm )
      set( RUN_SDC_M sdc_m__omp.slurm )
  elseif (OpenACC)
      set( RUN_SDC_K sdc_k__acc.slurm )
      set( RUN_SDC_M sdc_m__acc.slurm )
  else ()
      set( RUN_SDC_K sdc_k.slurm )
      set( RUN_SDC_M sdc_m.slurm )
  endif ()

  configure_file( ${RUN_SDC_K} sdc_k__slurm.sh )
  configure_file( ${RUN_SDC_M} sdc_m__slurm.sh )

endif (SLURM_FOUND)
