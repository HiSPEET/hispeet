file( COPY isp_flow__sdc_test.prm
           poiseuille.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (Slurm_FOUND)

  if (OpenMP)
  else ()
      configure_file( isp_flow__sdc_test.slurm isp_flow__sdc_test.slurm )
  endif ()

endif (Slurm_FOUND)
