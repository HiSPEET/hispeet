#===============================================================================
# ml_test_functionality

add_executable( ml_test_functionality
                ml_test_functionality.f )

target_link_libraries(ml_test_functionality HiBase)

file( COPY 
           ml_test_functionality.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR}
    )

#===============================================================================
# ml_test_elliptic

add_executable( ml_test_elliptic
                ml_test_elliptic.f )

target_link_libraries(ml_test_elliptic HiBase)

file( COPY 
           ml_test_elliptic.prm
           ml_test_elliptic_1l.prm
           ml_test_elliptic_global-p.prm
           ml_test_elliptic_global-h.prm
           ml_test_elliptic_local-h.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR}
    )
