add_executable(stencil_test_sol EXCLUDE_FROM_ALL stencil_test_sol.f)
#target_include_directories(stencil_test PUBLIC ~/HiSPEET/app/exercises/ex_5)
target_link_libraries(stencil_test_sol HiBase)
#target_include_directories(stencil_test PUBLIC $<TARGET_PROPERTY:Stencil_Baseline,Fortran_MODULE_DIRECTORY>)
#include(stencil_baseline)

add_executable(matrix_matrix_test_sol EXCLUDE_FROM_ALL matrix_matrix_test_sol.f)
target_link_libraries(matrix_matrix_test_sol HiBase)

add_executable(potpourri_a1_sol EXCLUDE_FROM_ALL potpourri_a1_sol.f)
target_link_libraries(potpourri_a1_sol HiBase)

add_executable(potpourri_a2_sol EXCLUDE_FROM_ALL potpourri_a2_sol.f)
target_link_libraries(potpourri_a2_sol HiBase)

add_executable(potpourri_a3_sol EXCLUDE_FROM_ALL potpourri_a3_sol.f)
target_link_libraries(potpourri_a3_sol HiBase)

add_executable(stencil_test EXCLUDE_FROM_ALL stencil_test.f)
#target_include_directories(stencil_test PUBLIC ~/HiSPEET/app/exercises/ex_5)
target_link_libraries(stencil_test HiBase)
#target_include_directories(stencil_test PUBLIC $<TARGET_PROPERTY:Stencil_Baseline,Fortran_MODULE_DIRECTORY>)
#include(stencil_baseline)

add_executable(matrix_matrix_test EXCLUDE_FROM_ALL matrix_matrix_test.f)
target_link_libraries(matrix_matrix_test HiBase)

add_executable(potpourri_a1 EXCLUDE_FROM_ALL potpourri_a1.f)
target_link_libraries(potpourri_a1 HiBase)

add_executable(potpourri_a2 EXCLUDE_FROM_ALL potpourri_a2.f)
target_link_libraries(potpourri_a2 HiBase)

add_executable(potpourri_a3 EXCLUDE_FROM_ALL potpourri_a3.f)
target_link_libraries(potpourri_a3 HiBase)

configure_file(README.md README.md COPYONLY)






