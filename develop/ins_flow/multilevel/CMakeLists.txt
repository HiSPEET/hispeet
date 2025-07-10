#===============================================================================
# ins_timeintegrator_3d_test

add_executable( ml_ins_solver_3d
                ml_ins_solver_3d.f
              )

target_link_libraries(ml_ins_solver_3d HiINSFlow)

#===============================================================================
# Test cases

add_subdirectory(hagen_poiseuille)
add_subdirectory(variable_viscosity)
