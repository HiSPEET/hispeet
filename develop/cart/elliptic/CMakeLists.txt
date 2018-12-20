#===============================================================================
# Elliptic test cases

# elliptic_test ................................................................

add_executable(elliptic_test elliptic_test.f elliptic_test_case.f)

target_link_libraries(elliptic_test HiCart)
configure_file(elliptic_test.prm elliptic_test.prm COPYONLY)

# elliptic_test_ci .............................................................

add_executable(elliptic_test elliptic_test__ip_ci.f elliptic_test_case.f)

target_link_libraries(elliptic_test__ip_ci HiCart)
configure_file(elliptic_test__ip_ci.prm elliptic_test__ip_ci.prm COPYONLY)

