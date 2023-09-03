add_executable(conservation_law
    conservation_law.f
    
    cl__operator__1d.f
    cl__problem__1d.f
    cl__problem__convection_diffusion__1d.f
    cl__problem__convection_diffusion__wave_package__1d.f
    cl__problem__burgers__1d.f
    cl__problem__burgers__moving_front__1d.f
    cl__problem__burgers__wave_package__1d.f
    
#   cl__problem__cns__1d               # NOT IMPLEMENTED YET!
#   cl__problem__cns__shock_tube__1d   # NOT IMPLEMENTED YET!

    cl__time_integrator__1d.f
    cl__time_integrator__euler__1d.f
    cl__time_integrator__isd1__1d.f
    cl__time_integrator__isd2__1d.f
    cl__time_integrator__rk__1d.f
   
    cl__sdc__method__1d.f
#   cl__sdc__method__euler__1d.f              # NOT IMPLEMENTED YET!
    cl__sdc__method__isd1__1d.f   

    cl__mlsdc__level__1d.f
    cl__mlsdc__level__1d__mp_interpolate_cf.f
    cl__mlsdc__level__1d__mp_interpolate_fc.f

    )

target_link_libraries(conservation_law HiBase HiPhysics)

#===============================================================================
# Configure files

#configure_file( time_integration.prm 
#                time_integration.prm COPYONLY)

configure_file( cl__problem__convection_diffusion__wave_package__1d.prm
                cl__problem__convection_diffusion__wave_package__1d.prm 
                COPYONLY )

configure_file( cl__problem__burgers__wave_package__1d.prm
                cl__problem__burgers__wave_package__1d.prm 
                COPYONLY )

configure_file( cl__problem__burgers__moving_front__1d.prm
                cl__problem__burgers__moving_front__1d.prm 
                COPYONLY )
