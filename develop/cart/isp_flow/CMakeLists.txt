#===============================================================================
# Executables

# isp_sdc__flow_test ...........................................................

add_executable(isp_flow__sdc_test isp_flow__sdc_test.F)

target_link_libraries(isp_flow__sdc_test HiCart)

# isp_sdc__compare_data ........................................................

add_executable(isp_flow__compare_data isp_flow__compare_data.F)

target_link_libraries(isp_flow__compare_data HiCart)

#===============================================================================
# Test cases

add_subdirectory(channel)
add_subdirectory(no_flow)
add_subdirectory(poiseuille)
add_subdirectory(stokes_dkm)
add_subdirectory(stokes_gms)
add_subdirectory(transition_tg)
add_subdirectory(vortex_hw)
add_subdirectory(vortex_tg)
add_subdirectory(variable_viscosity)
