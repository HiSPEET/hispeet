#===============================================================================
# ml_test_functionality

add_executable( ml_test_functionality
                ml_test_functionality.f )

target_link_libraries(ml_test_functionality HiBase)

configure_file( ml_test_functionality.prm
                ml_test_functionality.prm COPYONLY )

#===============================================================================
# ml_test_elliptic

add_executable( ml_test_elliptic
                ml_test_elliptic.f )

target_link_libraries(ml_test_elliptic HiBase)

configure_file( ml_test_elliptic.prm
                ml_test_elliptic.prm COPYONLY )
