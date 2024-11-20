#===============================================================================
# Incompressible Navier-Stokes flow library

add_library(HiINSFlow

            ins__operator__3d.f
            ins__operator__3d__mp_applydiffusionoperator_c.f
            ins__operator__3d__mp_applydiffusionoperator_v.f
            ins__operator__3d__mp_applyessentialbc.f
            ins__operator__3d__mp_applynaturalbc.f
            ins__operator__3d__mp_diffusionsolver.f
            ins__operator__3d__mp_getbackflowpenalty.f
            ins__operator__3d__mp_getdiffusionresidual_c.f
            ins__operator__3d__mp_getdiffusionresidual_v.f
            ins__operator__3d__mp_getdiffusionterm_c.F
            ins__operator__3d__mp_getdiffusionterm_v.F
            ins__operator__3d__mp_getstokesresidual.f
            ins__operator__3d__mp_getviscousboundarystress_c.f
            ins__operator__3d__mp_getviscousboundarystress_v.f
            ins__operator__3d__mp_pressuresolver.f

            ins__problem__3d.f
            ins__problem__channel__3d.f
            ins__problem__cylinder_2d__3d.f
            ins__problem__hagen_poiseuille__3d.f
            ins__problem__poiseuille__3d.f
            ins__problem__stokes_dkm__3d.f
            ins__problem__stokes_linke__3d.f
            ins__problem__vortex_tg__3d.f
            ins__problem__variable_viscosity__3d.f
            ins__problem__test_suite__3d.f

            ins__time_scales__3d.f

            ins__time_integrator__3d.f
            ins__time_integrator__3d__mp_fgmres_step.f
            ins__time_integrator__3d__mp_projectionstep.f
            ins__time_integrator__euler__3d.f
            ins__time_integrator__bdf2__3d.f
            ins__time_integrator__runge_kutta__3d.f

            ins__flow_characteristics__3d.f
            )


set_property(TARGET HiINSFlow
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(HiINSFlow PUBLIC HiBase)
