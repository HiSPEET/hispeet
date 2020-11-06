if (Slurm_FOUND)

  file( COPY isp_flow__sdc_test.tmpl
             variable_viscosity.prm
             convergence_dx.sh
        DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

  file(GLOB Slurm_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.slurm")
  foreach(FILE ${Slurm_FILES})
      configure_file( ${FILE} ${FILE} @ONLY )
  endforeach ()

endif (Slurm_FOUND)
