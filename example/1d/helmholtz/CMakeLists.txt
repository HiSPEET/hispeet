#===============================================================================
# Executables

# cg_helmholtz_1d ..............................................................

add_executable(cg_helmholtz_1d cg_helmholtz_1d.f)

target_link_libraries(cg_helmholtz_1d HiBase)

configure_file(cg_helmholtz_1d.prm cg_helmholtz_1d.prm COPYONLY)

# cg_helmholtz_1d_acc ..........................................................

add_executable(cg_helmholtz_1d_acc cg_helmholtz_1d_acc.f)

target_link_libraries(cg_helmholtz_1d_acc HiBase)

configure_file(cg_helmholtz_1d_acc.prm cg_helmholtz_1d_acc.prm COPYONLY)
