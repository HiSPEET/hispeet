#===============================================================================
# dg_elliptic_3d_test

add_executable( dg_elliptic_3d_test
                dg_elliptic_3d_test.f
                elliptic_problem.f
                elliptic_problem__knotty.f
                elliptic_problem__simple_1d.f
                elliptic_problem__simple_2d.f
                elliptic_problem__simple_3d.f
              )

target_link_libraries(dg_elliptic_3d_test HiBase)

configure_file( dg_elliptic_3d_test.prm
                dg_elliptic_3d_test.prm COPYONLY )

configure_file( dg_elliptic_3d_test.slurm
                dg_elliptic_3d_test.slurm COPYONLY )

configure_file( dg_elliptic_3d_test.tmpl
                dg_elliptic_3d_test.tmpl COPYONLY )

