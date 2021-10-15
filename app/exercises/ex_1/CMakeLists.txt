add_executable(hello EXCLUDE_FROM_ALL hello.f)

add_executable(calculate_pi EXCLUDE_FROM_ALL calculate_pi.f)
target_link_libraries(calculate_pi HiBase)

add_executable(jacobi_method EXCLUDE_FROM_ALL jacobi_method.f)
target_link_libraries(jacobi_method HiBase)
configure_file("jacobi_method.py" "jacobi_method.py" COPYONLY)