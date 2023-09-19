add_executable(convection_diffusion
               convection_diffusion.f
               cd__imex_bdf2__1d.f
               cd__imex_bdf3__1d.f
               cd__imex_euler__1d.f
               cd__imex_rk__1d.f
               cd__sdc_euler__1d.f
               cd__sdc_rk__1d.f
               cd__utils__1d.f
               )

target_link_libraries(convection_diffusion HiBase HiPhysics)

#===============================================================================
# Configure files

configure_file(convection_diffusion.prm  convection_diffusion.prm  COPYONLY)
configure_file(convergence_dt.py         convergence_dt.py         COPYONLY)
configure_file(convergence_dt.tmpl       convergence_dt.tmpl       COPYONLY)
configure_file(convergence_dt.sh         convergence_dt.sh         COPYONLY)

