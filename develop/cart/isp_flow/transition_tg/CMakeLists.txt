#===============================================================================
# Executables

add_executable( transition_tg_post transition_tg_post.F )
target_link_libraries( transition_tg_post HiCart )

#===============================================================================
# Files

file( COPY isp_flow__sdc_test.prm
           transition_tg.prm
           transition_tg_post.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (Slurm_FOUND)

    file( COPY test-cases/
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

endif (Slurm_FOUND)

