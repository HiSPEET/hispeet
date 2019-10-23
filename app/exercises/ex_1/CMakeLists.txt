add_executable(hello EXCLUDE_FROM_ALL hello.f)

add_executable(calculate_pi EXCLUDE_FROM_ALL calculate_pi.f)
target_link_libraries(calculate_pi HiBase)

add_executable(jacobi_iteration EXCLUDE_FROM_ALL jacobi_iteration.f)
target_link_libraries(jacobi_iteration HiBase)

configure_file(plot.py plot.py COPYONLY)

