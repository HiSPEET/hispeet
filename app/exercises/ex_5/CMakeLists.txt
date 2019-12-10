add_executable(stencil_test EXCLUDE_FROM_ALL stencil_test.f)
target_link_libraries(stencil_test HiBase)

add_executable(matrix_matrix_test EXCLUDE_FROM_ALL matrix_matrix_test.f)
target_link_libraries(matrix_matrix_test HiBase)

add_executable(potpourri_a1 EXCLUDE_FROM_ALL potpourri_a1.f)
target_link_libraries(potpourri_a1 HiBase)

add_executable(potpourri_a2 EXCLUDE_FROM_ALL potpourri_a2.f)
target_link_libraries(potpourri_a2 HiBase)

add_executable(potpourri_a3 EXCLUDE_FROM_ALL potpourri_a3.f)
target_link_libraries(potpourri_a3 HiBase)

configure_file(README.md README.md COPYONLY)






