#===============================================================================
# ins_timeintegrator_3d_test

add_executable( ins_timeintegrator_3d_test
                ins_timeintegrator_3d_test.f
              )

target_link_libraries(ins_timeintegrator_3d_test HiINSFlow)

#===============================================================================
# Test cases

add_subdirectory(no_flow)
add_subdirectory(stokes_linke)
add_subdirectory(vortex_tg)
add_subdirectory(variable_viscosity)

#===============================================================================
# Test cases - to be revised

file( COPY
          ins_timeintegrator_3d_test.prm
          channel
          cylinder_2d
          hagen_poiseuille
          poiseuille
          stokes_dkm
      DESTINATION
          ${CMAKE_CURRENT_BINARY_DIR} )

 configure_file( ins_timeintegrator_3d_test.slurm
                 ins_timeintegrator_3d_test.slurm COPYONLY )
