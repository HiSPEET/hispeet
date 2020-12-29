#===============================================================================
# mesh3d_explore

add_executable( mesh3d_explore
                mesh3d_explore.f )

target_link_libraries(mesh3d_explore HiBase)

#===============================================================================
# mesh3d_import_exodus

if (Exodus_FOUND)

    add_executable( mesh3d_import_exodus
                    mesh3d_import_exodus.f )

    target_link_libraries(mesh3d_import_exodus HiBase)

    #configure_file( mesh3d_import_exodus.prm
    #                mesh3d_import_exodus.prm COPYONLY )

    file( COPY exodus DESTINATION . )

endif ()
