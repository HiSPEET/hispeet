#===============================================================================
# Executables

# isp_sdc_flow_test ............................................................

add_executable(isp_flow__sdc_test isp_flow__sdc_test.f)

target_link_libraries(isp_flow__sdc_test HiCart)

#===============================================================================
# Test cases

add_subdirectory(channel)
add_subdirectory(poiseuille)
add_subdirectory(stokes_dkm)
add_subdirectory(stokes_gms)
add_subdirectory(transition_tg)
add_subdirectory(vortex_hw)
add_subdirectory(vortex_tg)
add_subdirectory(variable_viscosity)
