#===============================================================================
# mesh3d_explore

add_executable( mesh3d_explore
                mesh3d_explore.f )

target_link_libraries(mesh3d_explore HiBase)

configure_file( mesh3d_explore.prm
                mesh3d_explore.prm COPYONLY )

#===============================================================================
# mesh3d_import_generic

add_executable( mesh3d_import_generic
                mesh3d_import_generic.f )

target_link_libraries(mesh3d_import_generic HiBase)

configure_file( mesh3d_import_generic.prm
                mesh3d_import_generic.prm COPYONLY )

#===============================================================================
# mesh3d_import_gmsh

add_executable( mesh3d_import_gmsh
                mesh3d_import_gmsh.f )

target_link_libraries(mesh3d_import_gmsh HiBase)

configure_file( mesh3d_import_gmsh.prm
                mesh3d_import_gmsh.prm COPYONLY )

file( COPY gmsh DESTINATION . )

#===============================================================================
# mesh3d_test_orientation

add_executable( mesh3d_test_orientation
                mesh3d_test_orientation.f )

target_link_libraries(mesh3d_test_orientation HiBase)

configure_file( mesh3d_test_orientation.prm
                mesh3d_test_orientation.prm COPYONLY )

#===============================================================================
# mesh3d_adapt

add_executable( mesh3d_adapt
                mesh3d_adapt.f 
                smiling_face.f )

target_link_libraries(mesh3d_adapt HiBase)

configure_file( mesh3d_adapt.prm
                mesh3d_adapt.prm COPYONLY )
