file( COPY isp_flow__sdc_test.prm
           transition_tg.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (SLURM_FOUND)

    file( COPY test-cases/
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

endif (SLURM_FOUND)

