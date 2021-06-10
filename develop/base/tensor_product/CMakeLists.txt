#===============================================================================
# TPO AAA

add_executable( validate__tpo__aaa_3d
                validate__tpo__aaa_3d.f )

target_link_libraries(validate__tpo__aaa_3d HiBase)

configure_file( validate__tpo__aaa_3d.prm
                validate__tpo__aaa_3d.prm COPYONLY )

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
# TPO RotRot

add_executable( validate__tpo__rotrot_3d_r
                validate__tpo__rotrot_3d_r.f )

target_link_libraries(validate__tpo__rotrot_3d_r HiBase)

configure_file( validate__tpo__rotrot_3d_r.prm
                validate__tpo__rotrot_3d_r.prm COPYONLY )

#===============================================================================
# TPO Elliptic CI

add_executable( validate__tpo__elliptic_3d_rlci
                validate__tpo__elliptic_3d_rlci.f )

target_link_libraries(validate__tpo__elliptic_3d_rlci HiBase)

configure_file( validate__tpo__elliptic_3d_rlci.prm
                validate__tpo__elliptic_3d_rlci.prm COPYONLY )

#===============================================================================
# TPO Elliptic VI

add_executable( validate__tpo__elliptic_3d_rlvi
                validate__tpo__elliptic_3d_rlvi.f )

target_link_libraries(validate__tpo__elliptic_3d_rlvi HiBase)

configure_file( validate__tpo__elliptic_3d_rlvi.prm
                validate__tpo__elliptic_3d_rlvi.prm COPYONLY )

#===============================================================================
# TPO Spectral CI

add_executable( validate__tpo__spectral_3d_ci
                validate__tpo__spectral_3d_ci.f )

target_link_libraries(validate__tpo__spectral_3d_ci HiBase)

configure_file( validate__tpo__spectral_3d_ci.prm
                validate__tpo__spectral_3d_ci.prm COPYONLY )

#===============================================================================
# TPO Spectral CA

add_executable( validate__tpo__spectral_3d_ca
                validate__tpo__spectral_3d_ca.f )

target_link_libraries(validate__tpo__spectral_3d_ca HiBase)

configure_file( validate__tpo__spectral_3d_ca.prm
                validate__tpo__spectral_3d_ca.prm COPYONLY )

#===============================================================================
# TPO Schwarz CI

add_executable( validate__tpo__schwarz_3d_ci_rdp
                validate__tpo__schwarz_3d_ci_rdp.f )

target_link_libraries(validate__tpo__schwarz_3d_ci_rdp HiBase)

add_executable( validate__tpo__schwarz_3d_ci_rsp
                validate__tpo__schwarz_3d_ci_rsp.f )

target_link_libraries(validate__tpo__schwarz_3d_ci_rsp HiBase)

configure_file( validate__tpo__schwarz_3d_ci.prm
                validate__tpo__schwarz_3d_ci.prm COPYONLY )

#===============================================================================
# TPO Schwarz CA

add_executable( validate__tpo__schwarz_3d_ca
                validate__tpo__schwarz_3d_ca.f )

target_link_libraries(validate__tpo__schwarz_3d_ca HiBase)

configure_file( validate__tpo__schwarz_3d_ca.prm
                validate__tpo__schwarz_3d_ca.prm COPYONLY )

#===============================================================================
# Scripts

file( COPY scripts/validate__tpo.sh
           scripts/validate__tpo__aaa_3d.tmpl
           scripts/validate__tpo__div_3d_r.tmpl
           scripts/validate__tpo__elliptic_3d_rlci.tmpl
           scripts/validate__tpo__elliptic_3d_rlvi.tmpl
           scripts/validate__tpo__grad_3d_r.tmpl
           scripts/validate__tpo__rot_3d_r.tmpl
           scripts/validate__tpo__rotrot_3d_r.tmpl
           scripts/validate__tpo__schwarz_3d_ca.tmpl
           scripts/validate__tpo__schwarz_3d_ci.tmpl
           scripts/validate__tpo__spectral_3d_ca.tmpl
           scripts/validate__tpo__spectral_3d_ci.tmpl
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)

  if (OpenMP)
      set(RUN_SCRIPT scripts/validate__tpo_omp.slurm)
  elseif (OpenACC)
      set(RUN_SCRIPT scripts/validate__tpo_acc.slurm)
  else ()
      set(RUN_SCRIPT scripts/validate__tpo_seq.slurm)
  endif ()

  configure_file(${RUN_SCRIPT} validate__tpo.slurm)

endif (Slurm_FOUND)
