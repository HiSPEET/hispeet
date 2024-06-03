add_executable(kdtree_test kdtree_test.f)
target_link_libraries(kdtree_test KdTree)

add_executable(operatorspectrum operatorspectrum.f)
target_link_libraries(operatorspectrum HiBase)

add_executable(ip_spectrum ip_spectrum.f)
target_link_libraries(ip_spectrum HiBase)
configure_file(ip-input.prm ip-input.prm COPYONLY)

add_executable(radau radau.f)
target_link_libraries(radau HiBase)

add_executable(prep prep.F)
