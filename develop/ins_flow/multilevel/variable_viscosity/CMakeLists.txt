file( COPY test-diffusion-1g.prm
           test-diffusion-3g.prm
           test-diffusion-hp.prm
           test-projection-1g.prm
           test-projection-3g.prm
           test-projection-hp.prm
           test-solver-1g.prm
           test-solver-3g.prm
           test-solver-hp.prm
           variable_viscosity.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )
