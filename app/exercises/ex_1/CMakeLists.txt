add_executable(hello EXCLUDE_FROM_ALL hello.f)

add_executable(calculate_pi EXCLUDE_FROM_ALL calculate_pi.f)
target_link_libraries(calculate_pi HiBase)

add_executable(jacobi_iteration_sol EXCLUDE_FROM_ALL jacobi_iteration_sol.f)
target_link_libraries(jacobi_iteration_sol HiBase)

configure_file(plot_sol.py plot_sol.py COPYONLY)

configure_file(README.md README.md COPYONLY)
