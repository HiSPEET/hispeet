#===============================================================================
# ins_integrator_3d_test

add_executable( ins_integrator_3d_test
                ins_integrator_3d_test.f
              )

target_link_libraries(ins_integrator_3d_test HiINSFlow)

#===============================================================================
# Test cases

add_subdirectory(no_flow)
add_subdirectory(channel)
add_subdirectory(cylinder_2d)
add_subdirectory(hagen_poiseuille)
add_subdirectory(poiseuille)
add_subdirectory(stokes_dkm)
add_subdirectory(stokes_linke)
add_subdirectory(vortex_tg)
add_subdirectory(variable_viscosity)
