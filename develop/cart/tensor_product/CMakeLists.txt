#===============================================================================
# TPO Grad

add_executable( validate__cart__tpo_grad
                validate__cart__tpo_grad.f )

target_link_libraries(validate__cart__tpo_grad HiCart)

configure_file( validate__cart__tpo_grad.prm
                validate__cart__tpo_grad.prm COPYONLY )

#===============================================================================
# TPO Div

add_executable( validate__cart__tpo_div
                validate__cart__tpo_div.f )

target_link_libraries(validate__cart__tpo_div HiCart)

configure_file( validate__cart__tpo_div.prm
                validate__cart__tpo_div.prm COPYONLY )

#===============================================================================
# TPO Rot

add_executable( validate__cart__tpo_rot
                validate__cart__tpo_rot.f )

target_link_libraries(validate__cart__tpo_rot HiCart)

configure_file( validate__cart__tpo_rot.prm
                validate__cart__tpo_rot.prm COPYONLY )

#===============================================================================
# TPO Rot Rot

add_executable( validate__cart__tpo_rotrot
                validate__cart__tpo_rotrot.f )

target_link_libraries(validate__cart__tpo_rotrot HiCart)

configure_file( validate__cart__tpo_rotrot.prm
                validate__cart__tpo_rotrot.prm COPYONLY )

#===============================================================================
# TPO Elliptic CI

add_executable( validate__cart__tpo_elliptic_ci
                validate__cart__tpo_elliptic_ci.f )

target_link_libraries(validate__cart__tpo_elliptic_ci HiCart)

configure_file( validate__cart__tpo_elliptic_ci.prm
                validate__cart__tpo_elliptic_ci.prm COPYONLY )

#===============================================================================
# TPO Elliptic 3d RLCI

add_executable( validate__cart__tpo_elliptic_3d_rlci
                validate__cart__tpo_elliptic_3d_rlci.f )

target_link_libraries(validate__cart__tpo_elliptic_3d_rlci HiCart)

configure_file( validate__cart__tpo_elliptic_ci.prm
                validate__cart__tpo_elliptic_ci.prm COPYONLY )

#===============================================================================
# TPO Elliptic VI

add_executable( validate__cart__tpo_elliptic_vi
                validate__cart__tpo_elliptic_vi.f )

target_link_libraries(validate__cart__tpo_elliptic_vi HiCart)

configure_file( validate__cart__tpo_elliptic_vi.prm
                validate__cart__tpo_elliptic_vi.prm COPYONLY )

#===============================================================================
# TPO Spectral

add_executable( validate__cart__tpo_spectral
                validate__cart__tpo_spectral.f )

target_link_libraries(validate__cart__tpo_spectral HiCart)

configure_file( validate__cart__tpo_spectral.prm
                validate__cart__tpo_spectral.prm COPYONLY )

#===============================================================================
# TPO Spectral Iso

add_executable( validate__cart__tpo_spectral_iso
                validate__cart__tpo_spectral_iso.f )

target_link_libraries(validate__cart__tpo_spectral_iso HiCart)

configure_file( validate__cart__tpo_spectral_iso.prm
                validate__cart__tpo_spectral_iso.prm COPYONLY )

#===============================================================================
# TPO Schwarz

add_executable( validate__cart__tpo_schwarz
                validate__cart__tpo_schwarz.f )

target_link_libraries(validate__cart__tpo_schwarz HiCart)

configure_file( validate__cart__tpo_schwarz.prm
                validate__cart__tpo_schwarz.prm COPYONLY )

#===============================================================================
# TPO Schwarz Iso

add_executable( validate__cart__tpo_schwarz_iso
                validate__cart__tpo_schwarz_iso.f )

target_link_libraries(validate__cart__tpo_schwarz_iso HiCart)

configure_file( validate__cart__tpo_schwarz_iso.prm
                validate__cart__tpo_schwarz_iso.prm COPYONLY )

#===============================================================================
# Scripts

file( COPY scripts/validate__cart__tpo.sh
           scripts/validate__cart__tpo_div.prm.tmpl
           scripts/validate__cart__tpo_elliptic_ci.prm.tmpl
           scripts/validate__cart__tpo_elliptic_vi.prm.tmpl
           scripts/validate__cart__tpo_grad.prm.tmpl
           scripts/validate__cart__tpo_rot.prm.tmpl
           scripts/validate__cart__tpo_rotrot.prm.tmpl
           scripts/validate__cart__tpo_spectral_iso.prm.tmpl
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (SLURM_FOUND)

  if (OpenMP)
      set(RUN_SCRIPT scripts/validate__cart__tpo__slurm_omp.sh)
  elseif (OpenACC)
      set(RUN_SCRIPT scripts/validate__cart__tpo__slurm_acc.sh)
  else ()
      set(RUN_SCRIPT scripts/validate__cart__tpo__slurm_seq.sh)
  endif ()

  configure_file(${RUN_SCRIPT} validate__cart__tpo__slurm.sh )

endif (SLURM_FOUND)
