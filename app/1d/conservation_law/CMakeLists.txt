#===============================================================================
# CL Library

add_library( HiConservationLaw_1D 
             cl__operator__1d.f

             cl__problem__1d.f
             cl__problem__1d__sm__moment_limiter.f

             cl__problem__convection_diffusion__1d.f
             cl__problem__convection_diffusion__wave_package__1d.f

             cl__problem__burgers__1d.f
             cl__problem__burgers__moving_front__1d.f
             cl__problem__burgers__sine_wave__1d.f
             cl__problem__burgers__wave_package__1d.f

             cl__problem__cns__1d.f
             cl__problem__cns__1d__sm_diffusion.f
             cl__problem__cns__acoustic_wave__1d.f
             cl__problem__cns__contact_layer__1d.f
             cl__problem__cns__shock_tube__1d.f
             cl__problem__cns__shu_osher__1d.f

             cl__time_integrator__1d.f
             cl__time_integrator__euler__1d.f
             cl__time_integrator__isd1__1d.f
             cl__time_integrator__isd2__1d.f
             cl__time_integrator__rk__1d.f

             cl__sdc__method__1d.f
             cl__sdc__method__euler__1d.f
             cl__sdc__method__isd1__1d.f

             cl__mlsdc__1d.f
             cl__mlsdc__level__1d.f
             cl__mlsdc__level__1d__mp_applycorrector.f
             cl__mlsdc__level__1d__mp_applypredictor.f
             cl__mlsdc__level__1d__mp_getresidual.f
             cl__mlsdc__level__1d__mp_interpolate_cf.f
             cl__mlsdc__level__1d__mp_project_fc.f
             cl__mlsdc__level__1d__mp_restrict_fc.f
             cl__mlsdc__variable__1d.f
             cl__mlsdc__v_cycle__1d.f

           )

target_link_libraries(HiConservationLaw_1D PUBLIC HiBase HiPhysics)

#===============================================================================
# Single-level program

add_executable(conservation_law_1l conservation_law_1l.f)

target_link_libraries(conservation_law_1l HiConservationLaw_1D)

#===============================================================================
# Multi-level program

add_executable(conservation_law_ml conservation_law_ml.f)

target_link_libraries(conservation_law_ml HiConservationLaw_1D)

#===============================================================================
# Input files

file( COPY
          burgers
          convection_diffusion
          cns
      DESTINATION
          ${CMAKE_CURRENT_BINARY_DIR} )
