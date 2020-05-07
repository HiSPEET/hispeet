file( COPY isp_flow__sdc_test.prm
           transition_tg.prm
           runge_kutta_parameters.prm
           convergence_dt.py
           convergence_dt.sh
	   convergence_dt.slurm
           isp_flow__sdc_test.slurm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (SLURM_FOUND)

  if (OpenMP)
  else ()
      configure_file( isp_flow__sdc_test.slurm isp_flow__sdc_test.slurm )
  endif ()

endif (SLURM_FOUND)
