file( COPY README.md
           ins_timeintegrator_3d_test_1g.prm
           ins_timeintegrator_3d_test_hp-mg.prm
           ins_timeintegrator_3d_test_p-mg.prm
           variable_viscosity.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )
