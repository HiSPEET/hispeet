file( COPY README.md
           cylinder_2d.msh
           cylinder_2d.prm
           test_1g.prm
           test_mg_start.prm
           test_mg_restart.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )
