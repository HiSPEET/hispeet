file( COPY stokes_dkm.prm
           study.sh
           study.tmpl
           test.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)
    configure_file(study-mpi.slurm study.slurm)
endif (Slurm_FOUND)
