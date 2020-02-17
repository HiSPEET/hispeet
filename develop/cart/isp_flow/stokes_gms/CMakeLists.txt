file( COPY stokes_dkm.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

endif (SLURM_FOUND)
