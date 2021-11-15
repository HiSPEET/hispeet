#===============================================================================
# dg_diffusion3d_test

add_executable( dg_diffusion3d_test
                dg_diffusion3d_test.f
                create_cuboid_cartesian.f
                create_cuboid_diamonds.f
                create_cylinder.f
                create_annulus.f 
                elliptic_problem.f
                elliptic_problem__knotty.f
                elliptic_problem__simple_1d.f
                elliptic_problem__simple_2d.f
                elliptic_problem__simple_3d.f
              )

target_link_libraries(dg_diffusion3d_test HiBase)

configure_file( dg_diffusion3d_test.prm
                dg_diffusion3d_test.prm COPYONLY )

configure_file( dg_diffusion3d_test.slurm
                dg_diffusion3d_test.slurm COPYONLY )

