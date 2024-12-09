#===============================================================================
# ml_mesh_functionality

add_executable( ml_mesh_functionality
                ml_mesh_functionality.f )

target_link_libraries(ml_mesh_functionality HiBase)

file( COPY 
           ml_mesh_functionality.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR}
    )

#===============================================================================
# ml_spacetime_functionality

add_executable( ml_spacetime_functionality
                ml_spacetime_functionality.f )

target_link_libraries(ml_spacetime_functionality HiBase)

file( COPY 
           ml_spacetime_functionality_global_p.prm
           ml_spacetime_functionality_local_hp.prm
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
