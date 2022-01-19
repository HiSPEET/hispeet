#===============================================================================
# dg_diffusion3d_test

add_executable( ins_operator_3d_test
                ins_operator_3d_test.f
                create_cuboid_cartesian.f
                create_cuboid_diamonds.f
                create_cylinder.f
                create_annulus.f
              )

target_link_libraries(ins_operator_3d_test HiIncompressibleNS)

file( COPY
          ins_operator_3d_test.prm
          vortex_tg.prm
          variable_viscosity.prm
      DESTINATION
          ${CMAKE_CURRENT_BINARY_DIR} )

#configure_file( ins_operator_3d_test.slurm
#                ins_operator_3d_test.slurm COPYONLY )
#
#configure_file( ins_operator_3d_test.tmpl
#                ins_operator_3d_test.tmpl COPYONLY )
