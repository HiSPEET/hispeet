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

add_executable( validate__tpo__div_3d_d
                validate__tpo__div_3d_d.f )

target_link_libraries(validate__tpo__div_3d_d HiBase)

configure_file( validate__tpo__div_3d_d.prm
                validate__tpo__div_3d_d.prm COPYONLY )

#===============================================================================
# TPO Grad

add_executable( validate__tpo__grad_3d_r
                validate__tpo__grad_3d_r.f )

target_link_libraries(validate__tpo__grad_3d_r HiBase)

configure_file( validate__tpo__grad_3d_r.prm
                validate__tpo__grad_3d_r.prm COPYONLY )

add_executable( validate__tpo__grad_3d_d
                validate__tpo__grad_3d_d.f )

target_link_libraries(validate__tpo__grad_3d_d HiBase)

configure_file( validate__tpo__grad_3d_d.prm
                validate__tpo__grad_3d_d.prm COPYONLY )

#===============================================================================
# TPO INS_Convection

add_executable( validate__tpo__ins_convection_3d_d
                validate__tpo__ins_convection_3d_d.f )

target_link_libraries(validate__tpo__ins_convection_3d_d HiBase HiIncNS)

configure_file( validate__tpo__ins_convection_3d_d.prm
                validate__tpo__ins_convection_3d_d.prm COPYONLY )


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
# TPO Diffusion RLCI

add_executable( validate__tpo__diffusion_3d_rlci
                validate__tpo__diffusion_3d_rlci.f )

target_link_libraries(validate__tpo__diffusion_3d_rlci HiBase)

configure_file( validate__tpo__diffusion_3d_rlci.prm
                validate__tpo__diffusion_3d_rlci.prm COPYONLY )

#===============================================================================
# TPO Diffusion RLVI

add_executable( validate__tpo__diffusion_3d_rlvi
                validate__tpo__diffusion_3d_rlvi.f )

target_link_libraries(validate__tpo__diffusion_3d_rlvi HiBase)

configure_file( validate__tpo__diffusion_3d_rlvi.prm
                validate__tpo__diffusion_3d_rlvi.prm COPYONLY )

#===============================================================================
# TPO Diffusion DLCI

add_executable( validate__tpo__diffusion_3d_dlci
                validate__tpo__diffusion_3d_dlci.f )

target_link_libraries(validate__tpo__diffusion_3d_dlci HiBase)

configure_file( validate__tpo__diffusion_3d_dlci.prm
                validate__tpo__diffusion_3d_dlci.prm COPYONLY )

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

file( COPY scripts/omp_auto.sh
           scripts/validate__tpo_r.sh
           scripts/validate__tpo_d.sh
           scripts/validate__tpo_d_conv.sh
           scripts/validate__tpo__aaa_3d.tmpl
           scripts/validate__tpo__div_3d_d.tmpl
           scripts/validate__tpo__div_3d_r.tmpl
           scripts/validate__tpo__diffusion_3d_dlci.tmpl
           scripts/validate__tpo__diffusion_3d_rlci.tmpl
           scripts/validate__tpo__diffusion_3d_rlvi.tmpl
           scripts/validate__tpo__grad_3d_d.tmpl
           scripts/validate__tpo__grad_3d_r.tmpl
           scripts/validate__tpo__ins_convection_3d_d.tmpl
           scripts/validate__tpo__rot_3d_r.tmpl
           scripts/validate__tpo__rotrot_3d_r.tmpl
           scripts/validate__tpo__schwarz_3d_ca.tmpl
           scripts/validate__tpo__schwarz_3d_ci.tmpl
           scripts/validate__tpo__spectral_3d_ca.tmpl
           scripts/validate__tpo__spectral_3d_ci.tmpl
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

if (Slurm_FOUND)

    set(RUN_SCRIPT_R scripts/validate__tpo_r.slurm)
    set(RUN_SCRIPT_D scripts/validate__tpo_d.slurm)

    configure_file(${RUN_SCRIPT_R} validate__tpo_r.slurm)
    configure_file(${RUN_SCRIPT_D} validate__tpo_d.slurm)

endif (Slurm_FOUND)
