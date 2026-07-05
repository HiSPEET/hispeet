#===============================================================================
# Incompressible Navier-Stokes flow library

add_library(HiINSFlow

            operator/ins__operator__3d.f
            operator/ins__operator__3d__mp_applydiffusionoperator_c.f
            operator/ins__operator__3d__mp_applydiffusionoperator_v.f
            operator/ins__operator__3d__mp_applyessentialbc.f
            operator/ins__operator__3d__mp_applynaturalbc.f
            operator/ins__operator__3d__mp_applystokesoperator.f
            operator/ins__operator__3d__mp_diffusionsolver.f
            operator/ins__operator__3d__mp_getbackflowpenalty.f
            operator/ins__operator__3d__mp_getdiffusionresidual_c.f
            operator/ins__operator__3d__mp_getdiffusionresidual_v.f
            operator/ins__operator__3d__mp_getdiffusionterm_c.F
            operator/ins__operator__3d__mp_getdiffusionterm_v.F
            operator/ins__operator__3d__mp_getpressureboundaryvalues.f
            operator/ins__operator__3d__mp_getstokesresidual.f
            operator/ins__operator__3d__mp_getvariableviscosity.f
            operator/ins__operator__3d__mp_getviscousboundarystress_c.f
            operator/ins__operator__3d__mp_getviscousboundarystress_v.f
            operator/ins__operator__3d__mp_pressuresolver.f
            operator/ins__operator__3d__mp_stokesfgmres.f
            operator/ins__operator__3d__mp_stokesprojection.f

            problem/ins__problem__3d.f
            problem/ins__problem__channel__3d.f
            problem/ins__problem__cylinder_2d__3d.f
            problem/ins__problem__hagen_poiseuille__3d.f
            problem/ins__problem__linear_cascade__3d.f
            problem/ins__problem__no_flow__3d.f
            problem/ins__problem__poiseuille__3d.f
            problem/ins__problem__stokes_dkm__3d.f
            problem/ins__problem__stokes_linke__3d.f
            problem/ins__problem__vortex_tg__3d.f
            problem/ins__problem__variable_viscosity__3d.f
            problem/ins__problem__test_suite__3d.f

            integrator/ins__integrator__3d.f
            integrator/ins__integrator__bdf2__prepstep__3d.f
            integrator/ins__integrator__bdf2__3d.f
            integrator/ins__integrator__euler__3d.f
            integrator/ins__integrator__runge_kutta__3d.f

            util/ins__flow_characteristics__3d.f
            util/ins__time_scales__3d.f
            
            multilevel/ml__ins__flow_characteristics.f
            multilevel/ml__ins__operator__3d.f
            multilevel/ml__ins__operator__3d__mp_diffusionstep.f
            multilevel/ml__ins__operator__3d__mp_projectionstep.f
            multilevel/ml__ins__operator__3d__mp_mg_stokes_start.f
            multilevel/ml__ins__operator__3d__mp_mg_stokes_cycle.f
            multilevel/ml__ins__integrator__3d.f
            multilevel/ml__ins__integrator__bdf2__3d.f
            multilevel/ml__ins__time_scales__3d.f

            )


set_property(TARGET HiINSFlow
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(HiINSFlow PUBLIC HiBase)
