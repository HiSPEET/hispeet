#===============================================================================
# test_projection

add_executable( test_projection
                test_projection.f )

target_link_libraries(test_projection HiBase)

configure_file( test_projection.prm 
                test_projection.prm COPYONLY )

configure_file( test_projection.py 
                test_projection.py COPYONLY )

#===============================================================================
# test_hp_coarsening

add_executable( test_hp_coarsening
                test_hp_coarsening.f )

target_link_libraries(test_hp_coarsening HiBase)

configure_file( test_hp_coarsening_mode1.prm 
                test_hp_coarsening_mode1.prm COPYONLY )

configure_file( test_hp_coarsening_mode2.prm 
                test_hp_coarsening_mode2.prm COPYONLY )

configure_file( test_hp_coarsening.py 
                test_hp_coarsening.py COPYONLY )
