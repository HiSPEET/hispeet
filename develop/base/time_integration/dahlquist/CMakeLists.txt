add_executable( dahlquist 
                dahlquist.f
                dq__time_integrator.f
                dq__time_integrator__euler.f
                dq__time_integrator__isd.f
                dq__time_integrator__rk.f
                dq__time_integrator__tr.f
                dq__sdc__method.f
                dq__sdc__method__euler.f
                dq__sdc__method__rk.f
              )
              
target_link_libraries(dahlquist HiBase)

configure_file(dahlquist.prm dahlquist.prm COPYONLY)

file( COPY 
           eval_stability.py
           eval_accuracy.py
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )
                      