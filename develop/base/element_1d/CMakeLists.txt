#===============================================================================
# test_l2_projection

add_executable( test_l2_projection
                test_l2_projection.f )

target_link_libraries(test_l2_projection HiBase)

configure_file( test_l2_projection.prm 
                test_l2_projection.prm COPYONLY )

#===============================================================================
# test_hp_coarsening

add_executable( test_hp_coarsening
                test_hp_coarsening.f )

target_link_libraries(test_hp_coarsening HiBase)

configure_file( test_hp_coarsening_mode1.prm 
                test_hp_coarsening_mode1.prm COPYONLY )

configure_file( test_hp_coarsening_mode2.prm 
                test_hp_coarsening_mode2.prm COPYONLY )
