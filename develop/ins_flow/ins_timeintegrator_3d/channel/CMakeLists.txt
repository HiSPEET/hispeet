file( COPY channel.prm
           test_1g_start.prm
           test_1g_restart.prm
           test_mg_start.prm
           test_mg_restart.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)
    configure_file(test_channel.slurm test_channel.slurm COPYONLY)
endif (Slurm_FOUND)
