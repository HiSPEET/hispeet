file( COPY transition_tg.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

#if (Slurm_FOUND)
#    configure_file(test_mg_small.slurm test_mg_small.slurm COPYONLY)
#    configure_file(test_mg_large.slurm test_mg_large.slurm COPYONLY)
#    configure_file(test_mg_huge.slurm  test_mg_huge.slurm  COPYONLY)
#endif (Slurm_FOUND)

