#===============================================================================
# ins_timeintegrator_3d_test

add_executable( ml_ins_solver_3d
                ml_ins_solver_3d.f
                ml__ins__time_averaging__3d.f
              )

target_link_libraries(ml_ins_solver_3d HiINSFlow)

#===============================================================================
# Test cases

add_subdirectory(cylinder_2d)
add_subdirectory(hagen_poiseuille)
add_subdirectory(poiseuille)
add_subdirectory(variable_viscosity)
