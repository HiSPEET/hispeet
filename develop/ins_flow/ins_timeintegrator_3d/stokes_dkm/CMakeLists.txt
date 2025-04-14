file( COPY convergence_dx.sh
           stokes_dkm.prm
           test.prm
           test_convergence.tmpl
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)
    configure_file(convergence_dx-mpi.slurm convergence_dx.slurm)
endif (Slurm_FOUND)
