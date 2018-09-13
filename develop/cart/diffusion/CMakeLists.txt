#===============================================================================
# Executables

# elliptic_test ................................................................

add_executable(elliptic_test elliptic_test.f elliptic_test_case.f)

target_link_libraries(elliptic_test hibase hicart)

#===============================================================================
# Configure files

configure_file(elliptic_test.prm elliptic_test.prm COPYONLY)
