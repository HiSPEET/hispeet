file( COPY README.md
           linear_cascade.prm
           test_1g.prm
           test_mg_start.prm
           test_mg_restart.prm
           test_mg.tmpl
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)
    configure_file(test_mg.slurm test_mg.slurm COPYONLY)
endif (Slurm_FOUND)
