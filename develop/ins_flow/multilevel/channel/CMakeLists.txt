file(COPY channel.prm  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (Slurm_FOUND)
    configure_file(test_hp_re395.tmpl  test_hp_re395.tmpl  COPYONLY)
    configure_file(test_hp_re590.tmpl  test_hp_re590.tmpl  COPYONLY)
    configure_file(test_hp.slurm       test_hp.slurm       COPYONLY)
else ()
    configure_file(test_hp.prm   test_hp.prm   COPYONLY)
endif (Slurm_FOUND)
