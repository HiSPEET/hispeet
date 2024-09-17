#===============================================================================
# dg_elliptic_3d_test

add_executable( dg_elliptic_3d_test
                dg_elliptic_3d_test.f )

target_link_libraries(dg_elliptic_3d_test HiBase)

configure_file( dg_elliptic_3d_test_simple.prm
                dg_elliptic_3d_test_simple.prm COPYONLY )

configure_file( dg_elliptic_3d_test_knotty.prm
                dg_elliptic_3d_test_knotty.prm COPYONLY )

configure_file( dg_elliptic_3d_test_tgv_pressure.prm
                dg_elliptic_3d_test_tgv_pressure.prm COPYONLY )

configure_file( dg_elliptic_3d_test.slurm
                dg_elliptic_3d_test.slurm COPYONLY )

configure_file( dg_elliptic_3d_test.tmpl
                dg_elliptic_3d_test.tmpl COPYONLY )

