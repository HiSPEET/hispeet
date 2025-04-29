file( COPY hagen_poiseuille.prm
           pipe.msh
           test_gmsh.prm
           test_inbuilt.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )
