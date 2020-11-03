file( COPY isp_flow__sdc_test.tmpl
           isp_flow__sdc_test.prm
           vortex_tg.prm
           sdc_k.sh
           sdc_k.py
           sdc_m.sh
           sdc_m.py
           convergence_dt.py
           convergence_dt.sh
	   convergence_dt.slurm
	   isp_flow__sdc_test.slurm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (Slurm_FOUND)

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

  configure_file( ${RUN_SDC_K} sdc_k.slurm )
  configure_file( ${RUN_SDC_M} sdc_m.slurm )

endif (Slurm_FOUND)
