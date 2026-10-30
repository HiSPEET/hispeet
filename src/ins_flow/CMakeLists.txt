#===============================================================================
# Incompressible Navier-Stokes flow library

add_library(HiINSFlow

            problem/ins__problem__3d.f
            problem/ins__problem__channel__3d.f
            problem/ins__problem__cylinder_2d__3d.f
            problem/ins__problem__hagen_poiseuille__3d.f
            problem/ins__problem__linear_cascade__3d.f
            problem/ins__problem__no_flow__3d.f
            problem/ins__problem__poiseuille__3d.f
            problem/ins__problem__stokes_dkm__3d.f
            problem/ins__problem__stokes_linke__3d.f
            problem/ins__problem__transition_tg__3d.f
            problem/ins__problem__vortex_tg__3d.f
            problem/ins__problem__variable_viscosity__3d.f
            problem/ins__problem__test_suite__3d.f

            sgs_model/ins__sgs__model__3d.f

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

            integrator/ins__integrator__3d.f
            integrator/ins__integrator__bdf__prepstep__3d.f
            integrator/ins__integrator__bdf__3d.f
            integrator/ins__integrator__euler__3d.f
            integrator/ins__integrator__runge_kutta__3d.f

            util/ins__boundary_fluxes__3d.f
            util/ins__flow_characteristics__3d.f
            util/ins__time_scales__3d.f
            
            multilevel/ml__ins__boundary_fluxes.f
            multilevel/ml__ins__flow_characteristics.f
            multilevel/ml__ins__operator__3d.f
            multilevel/ml__ins__diffusion__3d.f
            multilevel/ml__ins__diffusion__3d_mp_diffusionresidual.f
            multilevel/ml__ins__diffusion__3d_mp_diffusionstep.f
            multilevel/ml__ins__projection__3d.f
            multilevel/ml__ins__stokes__3d.f
            multilevel/ml__ins__stokes__3d__mp_stokescycle.f
            multilevel/ml__ins__stokes__3d__mp_stokesresidual.f
            multilevel/ml__ins__stokes__3d__mp_stokesstart.f
            multilevel/ml__ins__stokes__3d__mp_stokesstep.f
            multilevel/ml__ins__integrator__3d.f
            multilevel/ml__ins__integrator__bdf__3d.f
            multilevel/ml__ins__time_scales__3d.f

            )


set_property(TARGET HiINSFlow
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(HiINSFlow PUBLIC HiBase)

#===============================================================================
# User-defined problem
#
#   -D INS_USER_PROBLEM=<name>|<directory>
#
# selects the directory user/ins_problem/<name> or an absolute directory path.
# Its CMakeLists.txt calls hispeet_ins_user_problem(), which adds the sources
# of module INS__Problem__User__3D, and possibly user libraries, to HiINSFlow.
# Without user problem, a dummy is compiled. See user/ins_problem/README.md.

set(INS_USER_PROBLEM "" CACHE STRING
    "user-defined INS problem: name in user/ins_problem or directory")

#-------------------------------------------------------------------------------
# hispeet_ins_user_problem(SOURCES src... [LIBRARIES lib...])
#
#   SOURCES    user problem sources, relative to the calling directory
#   LIBRARIES  libraries used by the problem, e.g. imported targets from
#              find_package()

function(hispeet_ins_user_problem)
    cmake_parse_arguments(ARG "" "" "SOURCES;LIBRARIES" ${ARGN})
    if (ARG_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "hispeet_ins_user_problem: invalid arguments "
                            "'${ARG_UNPARSED_ARGUMENTS}'")
    endif ()
    if (NOT ARG_SOURCES)
        message(FATAL_ERROR "hispeet_ins_user_problem: no SOURCES given")
    endif ()
    target_sources(HiINSFlow PRIVATE ${ARG_SOURCES})
    if (ARG_LIBRARIES)
        target_link_libraries(HiINSFlow PUBLIC ${ARG_LIBRARIES})
    endif ()
endfunction()

#-------------------------------------------------------------------------------
# selection

if (INS_USER_PROBLEM)
    if (IS_ABSOLUTE "${INS_USER_PROBLEM}")
        set(user_dir "${INS_USER_PROBLEM}")
    else ()
        set(user_dir
            "${PROJECT_SOURCE_DIR}/user/ins_problem/${INS_USER_PROBLEM}")
    endif ()
    if (NOT EXISTS "${user_dir}/CMakeLists.txt")
        message(FATAL_ERROR "INS_USER_PROBLEM: no CMakeLists.txt in ${user_dir}")
    endif ()
    message(STATUS "INS user problem: ${user_dir}")
    add_subdirectory("${user_dir}" "${CMAKE_CURRENT_BINARY_DIR}/user_problem")
else ()
    message(STATUS "INS user problem: none (dummy)")
    target_sources(HiINSFlow PRIVATE problem/ins__problem__user__3d.f)
endif ()
