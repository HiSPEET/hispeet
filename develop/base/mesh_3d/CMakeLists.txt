#===============================================================================
# mesh3d_adapt

add_executable( mesh3d_adapt
                mesh3d_adapt.f 
                qoi__distribution__3d.f
                qoi__point__3d.f
                qoi__sphere__3d.f
                qoi__smiley__3d.f       )

target_link_libraries(mesh3d_adapt HiBase)

configure_file( mesh3d_adapt_point.prm
                mesh3d_adapt_point.prm COPYONLY )

configure_file( mesh3d_adapt_sphere.prm
                mesh3d_adapt_sphere.prm COPYONLY )

configure_file( mesh3d_adapt_smiley.prm
                mesh3d_adapt_smiley.prm COPYONLY )

#===============================================================================
# mesh3d_explore

add_executable( mesh3d_explore
                mesh3d_explore.f )

target_link_libraries(mesh3d_explore HiBase)

configure_file( mesh3d_explore.prm
                mesh3d_explore.prm COPYONLY )

#===============================================================================
# mesh3d_hdf5

add_executable( mesh3d_hdf5
                mesh3d_hdf5.f 
                qoi__distribution__3d.f
                qoi__smiley__3d.f       )

target_link_libraries(mesh3d_hdf5 HiBase)

configure_file( mesh3d_hdf5.prm
                mesh3d_hdf5.prm COPYONLY )

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

#===============================================================================
# mesh3d_test_orientation

add_executable( mesh3d_test_orientation
                mesh3d_test_orientation.f )

target_link_libraries(mesh3d_test_orientation HiBase)

configure_file( mesh3d_test_orientation.prm
                mesh3d_test_orientation.prm COPYONLY )
