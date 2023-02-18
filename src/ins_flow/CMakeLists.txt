#===============================================================================
# Incompressible Navier-Stokes flow library

add_library(HiINSFlow
            ins__operator__3d.f
            ins__operator__3d__mp_applyvelocitybc.f
            ins__operator__3d__mp_applydiffusionoperator_c.f
            ins__operator__3d__mp_diffusionsolver.f
            ins__operator__3d__mp_getdiffusionterm_dc.F
            ins__operator__3d__mp_getdiffusionresidual_c.f
            ins__operator__3d__mp_pressuresolver.f
            ins__problem__3d.f
            ins__problem__vortex_tg__3d.f
            ins__problem__variable_viscosity__3d.f
            ins__time_integrator__3d.f
            ins__time_integrator__3d__mp_projectionstep.f
            ins__time_integrator__euler__3d.f
            ins__time_integrator__bdf2__3d.f
            ins__flow_characteristics__3d.f
            )


set_property(TARGET HiINSFlow
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(HiINSFlow PUBLIC HiBase)
