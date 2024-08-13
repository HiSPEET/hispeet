#===============================================================================
# all

file( COPY gmsh DESTINATION . )

#===============================================================================
# ml_test_functionality

add_executable( ml_test_functionality
                ml_test_functionality.f )

target_link_libraries(ml_test_functionality HiBase)

configure_file( ml_test_functionality.prm
                ml_test_functionality.prm COPYONLY )
