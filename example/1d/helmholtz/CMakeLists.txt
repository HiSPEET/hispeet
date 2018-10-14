#===============================================================================
# Executables

# cg_helmholtz_1d ..............................................................

add_executable(cg_helmholtz_1d cg_helmholtz_1d.f)

target_link_libraries(cg_helmholtz_1d HiBase)

#===============================================================================
# Configure files

configure_file(cg_helmholtz_1d.prm cg_helmholtz_1d.prm COPYONLY)
