#===============================================================================
# Incompressible Navier-Stokes library

add_library(HiIncompressibleNS 
            ins__operator__3d.f
            ins__operator__3d__mp_getdiffusionterm_dc.f
            ins__problem__3d.f
            ins__problem__3d__vortex_tg.f
            ins__problem__3d__variable_viscosity.f
            )

set_property(TARGET HiIncompressibleNS
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(HiIncompressibleNS PUBLIC HiBase)
