file( COPY channel.prm
           test_1g_start.prm
           test_1g_restart.prm
           test_mg.tmpl
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)
    configure_file(test_mg.slurm test_mg.slurm COPYONLY)
endif (Slurm_FOUND)
