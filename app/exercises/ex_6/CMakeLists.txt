configure_file(input input COPYONLY)

add_executable(test_interpolation EXCLUDE_FROM_ALL test_interpolation.f spectral_element_interpolation.f)
target_link_libraries(test_interpolation HiBase)

add_executable(potential_flow_solver EXCLUDE_FROM_ALL potential_flow_solver.f helmholtz_2d.f)
target_link_libraries(potential_flow_solver HiBase)
#configure_file(README.md README.md COPYONLY)

