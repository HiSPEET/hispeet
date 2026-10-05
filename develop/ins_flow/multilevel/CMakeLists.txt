#===============================================================================
# test_ml_ins_projection_3d

add_executable( test_ml_ins_projection_3d
                test_ml_ins_projection_3d.f
              )

target_link_libraries(test_ml_ins_projection_3d HiINSFlow)

#===============================================================================
# test_ml_ins_diffusion_3d

add_executable( test_ml_ins_diffusion_3d
                test_ml_ins_diffusion_3d.f
              )

target_link_libraries(test_ml_ins_diffusion_3d HiINSFlow)

#===============================================================================
# ml_ins_solver_3d

# add_executable( ml_ins_solver_3d
#                 ml_ins_solver_3d.f
#                 ml__ins__time_averaging__3d.f
#               )
# 
# target_link_libraries(ml_ins_solver_3d HiINSFlow)

#===============================================================================
# Test cases

add_subdirectory(channel)
add_subdirectory(cylinder_2d)
add_subdirectory(hagen_poiseuille)
add_subdirectory(poiseuille)
add_subdirectory(variable_viscosity)
