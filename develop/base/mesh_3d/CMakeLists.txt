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
# mesh3d_import_exodus

if (Exodus_FOUND)

    add_executable( mesh3d_import_exodus
                    mesh3d_import_exodus.f )

    target_link_libraries(mesh3d_import_exodus HiBase)

    file( COPY exodus DESTINATION . )

endif ()

#===============================================================================
# mesh3d_test_orientation

add_executable( mesh3d_test_orientation
                mesh3d_test_orientation.f )

target_link_libraries(mesh3d_test_orientation HiBase)

configure_file( mesh3d_test_orientation.prm
                mesh3d_test_orientation.prm COPYONLY )

