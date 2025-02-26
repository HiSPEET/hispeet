file( COPY channel.prm
           test_1g_start.prm
           test_1g_restart.prm
           test_mg_large.tmpl
           test_mg_small.tmpl
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)
    configure_file(test_mg_large.slurm test_mg_large.slurm COPYONLY)
    configure_file(test_mg_small.slurm test_mg_small.slurm COPYONLY)
endif (Slurm_FOUND)

