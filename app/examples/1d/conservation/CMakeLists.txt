add_executable(dg_conservation_1d
    dg_conservation_1d.f
    conservation_problem_1d.f
    scalar_problem_1d.f
    burgers_problem_1d.f
    burgers_problem_1d__breaking_wave.f
    )

target_link_libraries(dg_conservation_1d HiBase)

#===============================================================================
# Configure files

configure_file( dg_conservation_1d.prm 
                dg_conservation_1d.prm COPYONLY )

configure_file( burgers_problem_1d__breaking_wave.prm 
                burgers_problem_1d__breaking_wave.prm COPYONLY )

