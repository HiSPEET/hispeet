file( COPY ins_timeintegrator_3d_test_2d.prm
           ins_timeintegrator_3d_test_3d.prm
           no_flow_2d.prm
           no_flow_3d.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)
   configure_file(no_flow_2d.slurm no_flow_2d.slurm)
   configure_file(no_flow_3d.slurm no_flow_3d.slurm)
endif (Slurm_FOUND)
