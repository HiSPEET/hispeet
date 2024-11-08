#===============================================================================
# ins_timeintegrator_3d_test

add_executable( ins_timeintegrator_3d_test
                ins_timeintegrator_3d_test.f
              )

target_link_libraries(ins_timeintegrator_3d_test HiINSFlow)

#===============================================================================
# Test cases

add_subdirectory(vortex_tg)

#===============================================================================
# Test cases - to be revised

file( COPY
          ins_timeintegrator_3d_test.prm
          channel
          cylinder_2d
          hagen_poiseuille
          poiseuille
          stokes_dkm
          stokes_linke
          variable_viscosity
          #vortex_tg             # -- revised 
      DESTINATION
          ${CMAKE_CURRENT_BINARY_DIR} )

 configure_file( ins_timeintegrator_3d_test.slurm
                 ins_timeintegrator_3d_test.slurm COPYONLY )
