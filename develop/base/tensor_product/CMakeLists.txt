#===============================================================================
# TPO Div

add_executable( validate__tpo__div_3d_r
                validate__tpo__div_3d_r.f )

target_link_libraries(validate__tpo__div_3d_r HiBase)

configure_file( validate__tpo__div_3d_r.prm
                validate__tpo__div_3d_r.prm COPYONLY )

#===============================================================================
# TPO Grad

add_executable( validate__tpo__grad_3d_r
                validate__tpo__grad_3d_r.f )

target_link_libraries(validate__tpo__grad_3d_r HiBase)

configure_file( validate__tpo__grad_3d_r.prm
                validate__tpo__grad_3d_r.prm COPYONLY )

#===============================================================================
# TPO Rot

add_executable( validate__tpo__rot_3d_r
                validate__tpo__rot_3d_r.f )

target_link_libraries(validate__tpo__rot_3d_r HiBase)

configure_file( validate__tpo__rot_3d_r.prm
                validate__tpo__rot_3d_r.prm COPYONLY )

#===============================================================================
# TPO Elliptic CI

add_executable( validate__tpo__elliptic_3d_rlci
                validate__tpo__elliptic_3d_rlci.f )

target_link_libraries(validate__tpo__elliptic_3d_rlci HiBase)

configure_file( validate__tpo__elliptic_3d_rlci.prm
                validate__tpo__elliptic_3d_rlci.prm COPYONLY )

##===============================================================================
# Scripts
#
#file( COPY
#           scripts/validate__tpo__elliptic_3d_rlci.prm.tmpl
#      DESTINATION
#           ${CMAKE_CURRENT_BINARY_DIR} )
#
#if (SLURM_FOUND)
#
#  if (OpenMP)
#      set(RUN_SCRIPT scripts/validate__cart__tpo__slurm_omp.sh)
#  elseif (OpenACC)
#      set(RUN_SCRIPT scripts/validate__cart__tpo__slurm_acc.sh)
#  else ()
#      set(RUN_SCRIPT scripts/validate__cart__tpo__slurm_seq.sh)
#  endif ()
#
#  configure_file(${RUN_SCRIPT} validate__cart__tpo__slurm.sh )
#
#endif (SLURM_FOUND)
