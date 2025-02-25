#===============================================================================
# ins_timeintegrator_3d_test

add_executable( ins_timeintegrator_3d_test
                ins_timeintegrator_3d_test.f
              )

target_link_libraries(ins_timeintegrator_3d_test HiINSFlow)

#===============================================================================
# Test cases

add_subdirectory(no_flow)
add_subdirectory(channel)
add_subdirectory(hagen_poiseuille)
add_subdirectory(poiseuille)
add_subdirectory(stokes_linke)
add_subdirectory(vortex_tg)
add_subdirectory(variable_viscosity)

#===============================================================================
# Test cases - to be revised

file( COPY
          cylinder_2d
          stokes_dkm
      DESTINATION
          ${CMAKE_CURRENT_BINARY_DIR} )
