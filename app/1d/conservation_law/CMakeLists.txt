add_executable(conservation_law
    conservation_law.f
    
    cl__problem__1d.f
    cl__problem__scalar__1d.f
    cl__problem__scalar__burgers__1d.f
    cl__problem__scalar__burgers__breaking_wave__1d.f
    
#   cl__problem__system__1d                    # NOT IMPLEMENTED YET!
#   cl__problem__system__cns__1d               # NOT IMPLEMENTED YET!
#   cl__problem__system__cns__shock_tube__1d   # NOT IMPLEMENTED YET!

    cl__time_integrator__1d.f
    cl__time_integrator__tvdrk__1d.f
    cl__time_integrator__euler__1d.f
#   cl__time_integrator__isd__1d.f             # NOT IMPLEMENTED YET!
    cl__time_integrator__rk__1d.f
    cl__time_integrator__tr__1d.f
   
    cl__sdc__method__1d.f
    cl__sdc__method__euler__1d.f
#   cl__sdc__method__isd__1d.f                 # NOT IMPLEMENTED YET!
    cl__sdc__method__rk__1d.f
    cl__sdc__method__tvdrk__1d.f
    )

target_link_libraries(conservation_law HiBase)

#===============================================================================
# Configure files

configure_file( time_integration.prm 
                time_integration.prm COPYONLY)

configure_file( conservation_law.prm 
                conservation_law.prm COPYONLY )

configure_file( cl__problem__scalar__burgers__breaking_wave__1d.prm 
                cl__problem__scalar__burgers__breaking_wave__1d.prm COPYONLY )

