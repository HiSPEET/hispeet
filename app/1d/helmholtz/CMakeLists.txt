#===============================================================================
# Executables

# cg_helmholtz ..............................................................

add_executable( cg_helmholtz
                cg_helmholtz.f
                helmholtz_test_cases.f )

target_link_libraries(cg_helmholtz HiBase)

configure_file( cg_helmholtz.prm
                cg_helmholtz.prm COPYONLY )

# cg_helmholtz_acc ..........................................................

add_executable( cg_helmholtz_acc
                cg_helmholtz_acc.f
                helmholtz_test_cases.f )

target_link_libraries(cg_helmholtz_acc HiBase)

configure_file( cg_helmholtz_acc.prm
                cg_helmholtz_acc.prm COPYONLY )

# dg_helmholtz ..............................................................

add_executable( dg_helmholtz
                dg_helmholtz.f
                helmholtz_test_cases.f )

target_link_libraries(dg_helmholtz HiBase)

configure_file( dg_helmholtz.prm
                dg_helmholtz.prm COPYONLY )
