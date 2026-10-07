#===============================================================================
# Elliptic test cases

# elliptic_test_ci .............................................................

add_executable( elliptic_test__ip_ci
                elliptic_test__ip_ci.f
                elliptic_problem.f
                elliptic_problem__simple_1d.f
                elliptic_problem__simple_2d.f
                elliptic_problem__simple_3d.f
                elliptic_problem__knotty.f    )

target_link_libraries( elliptic_test__ip_ci HiCart )

configure_file( elliptic_test__ip_ci.prm
                elliptic_test__ip_ci.prm COPYONLY )

configure_file( elliptic_test__ip_ci.slurm
                elliptic_test__ip_ci.slurm )

# elliptic_test_vi .............................................................

add_executable( elliptic_test__ip_vi
                elliptic_test__ip_vi.f
                elliptic_problem.f
                elliptic_problem__simple_1d.f
                elliptic_problem__simple_2d.f
                elliptic_problem__simple_3d.f
                elliptic_problem__knotty.f    )

target_link_libraries( elliptic_test__ip_vi HiCart )

configure_file( elliptic_test__ip_vi.prm
                elliptic_test__ip_vi.prm COPYONLY )

#===============================================================================
# separate module directories: sources shared by several executables are
# compiled per target; a common module directory causes races in parallel builds

foreach (t elliptic_test__ip_ci elliptic_test__ip_vi)
    set_target_properties(${t} PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules_${t})
endforeach ()
