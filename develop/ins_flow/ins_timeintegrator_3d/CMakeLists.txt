#===============================================================================
# dg_diffusion3d_test

add_executable( ins_timeintegrator_3d_test
                ins_timeintegrator_3d_test.f
              )

target_link_libraries(ins_timeintegrator_3d_test HiINSFlow)

file( COPY
          ins_timeintegrator_3d_test.prm
          stokes_dkm
          stokes_linke
          variable_viscosity
          vortex_tg
      DESTINATION
          ${CMAKE_CURRENT_BINARY_DIR} )

 configure_file( ins_timeintegrator_3d_test.slurm
                 ins_timeintegrator_3d_test.slurm COPYONLY )
 
#configure_file( ins_operator_3d_test.tmpl
#                ins_operator_3d_test.tmpl COPYONLY )
