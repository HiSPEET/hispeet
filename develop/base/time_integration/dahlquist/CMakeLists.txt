#===============================================================================
# single-level

add_executable( dahlquist 
                dahlquist.f
                dq__time_integrator.f
                dq__sdc__method__picard.f
                dq__time_integrator__euler.f
                dq__time_integrator__isd.f
                dq__time_integrator__rk.f
                dq__time_integrator__tr.f
                dq__sdc__method.f
                dq__sdc__method__euler.f
                dq__sdc__method__diag.f
                dq__sdc__method__diag_sd.f
                dq__sdc__method__isd.f
                dq__sdc__method__lu.f
                dq__sdc__method__rk.f
              )
              
target_link_libraries(dahlquist HiBase)

configure_file(dahlquist.prm dahlquist.prm COPYONLY)

#===============================================================================
# multi-level

add_executable( dahlquist_mlsdc
                dahlquist_mlsdc.f
                dq__time_integrator.f
                dq__time_integrator__euler.f
                dq__time_integrator__isd.f
                dq__time_integrator__rk.f
                dq__time_integrator__tr.f
                dq__sdc__method.f
                dq__sdc__method__euler.f
                dq__sdc__method__isd.f
                dq__mlsdc__corrector__euler.f
                dq__mlsdc__corrector__euler_bepar.f
                dq__mlsdc__corrector__isd.f
                dq__mlsdc__1d.f
                dq__mlsdc__level__1d.f
                dq__mlsdc__level__1d__mp_applycorrector.f
                dq__mlsdc__level__1d__mp_applyoperator.f
                dq__mlsdc__level__1d__mp_applypredictor.f
                dq__mlsdc__level__1d__mp_getresidual.f
                dq__mlsdc__level__1d__mp_interpolate_cf.f
                dq__mlsdc__level__1d__mp_project_fc.f
                dq__mlsdc__level__1d__mp_restrict_fc.f
                dq__mlsdc__variable__1d.f
                dq__mlsdc__cascade__1d.f
                dq__mlsdc__fmg__1d.f
                dq__mlsdc__v_cycle__1d.f
              )

target_link_libraries(dahlquist_mlsdc HiBase)

configure_file(dahlquist_mlsdc.prm dahlquist_mlsdc.prm COPYONLY)

file( COPY 
           mlsdc_study 
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

#===============================================================================
# evaluation

file( COPY 
           eval_stability.py
           eval_accuracy.py
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )
