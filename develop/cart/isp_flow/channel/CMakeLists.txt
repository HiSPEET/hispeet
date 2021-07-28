#===============================================================================
# Executables

add_executable( channel_post channel_post.F )
target_link_libraries( channel_post HiCart )

add_executable( change_polynomial_order change_polynomial_order.f )
target_link_libraries( change_polynomial_order HiCart )

#===============================================================================
# Configuration files and scripts

file( COPY isp_flow__sdc_test.prm
           channel.prm
           channel_post.prm
           change_polynomial_order.prm
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )


if (Slurm_FOUND)

  if (OpenMP)
  else ()
      configure_file( isp_flow__sdc_test.slurm isp_flow__sdc_test.slurm )
  endif ()

endif (Slurm_FOUND)
