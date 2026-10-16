add_executable(stencil_test EXCLUDE_FROM_ALL stencil_test.f)
target_link_libraries(stencil_test HiBase)

add_executable(matrix_matrix_test EXCLUDE_FROM_ALL matrix_matrix_test.f)
target_link_libraries(matrix_matrix_test HiBase)

add_executable(potpourri_1 EXCLUDE_FROM_ALL potpourri_1.f potpourri_print_mlups.f)
target_link_libraries(potpourri_1 HiBase)

add_executable(potpourri_2 EXCLUDE_FROM_ALL potpourri_2.f potpourri_print_mlups.f)
target_link_libraries(potpourri_2 HiBase)

add_executable(potpourri_3 EXCLUDE_FROM_ALL potpourri_3.f potpourri_print_mlups.f)
target_link_libraries(potpourri_3 HiBase)

configure_file(README.md README.md COPYONLY)






