#===============================================================================
# ml_mesh_functionality

add_executable( ml_mesh_functionality
                ml_mesh_functionality.f )

target_link_libraries(ml_mesh_functionality HiBase)

file( COPY 
           ml_mesh_functionality_gmsh.prm
           ml_mesh_functionality_cuboid.prm
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
# ml_elliptic_test_static

add_executable( ml_elliptic_test_static
                ml_elliptic_test_static.f )

target_link_libraries(ml_elliptic_test_static HiBase)

file( COPY 
           ml_elliptic_test_static.prm
           ml_elliptic_test_static_1l.prm
           ml_elliptic_test_static_global-p.prm
           ml_elliptic_test_static_global-h.prm
           ml_elliptic_test_static_local-h.prm
           robustness
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR}
    )

#===============================================================================
# ml_elliptic_test_adaptive

add_executable( ml_elliptic_test_adaptive
                ml_elliptic_test_adaptive.f )

target_link_libraries(ml_elliptic_test_adaptive HiBase)

file( COPY 
           ml_elliptic_test_adaptive.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR}
    )

configure_file( ml_elliptic_test_adaptive.slurm
                ml_elliptic_test_adaptive.slurm COPYONLY )

