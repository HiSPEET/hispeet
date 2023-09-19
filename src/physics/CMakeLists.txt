#===============================================================================
# 1D

set(SCALAR_1D scalar_1d/harmonic_wave_package.f)

#===============================================================================
# isp_flow

file(GLOB ISP_FLOW isp_flow/*.f)

#===============================================================================
# physics library

set(PHYSICS_SRC ${SCALAR_1D}
                ${ISP_FLOW}
                )

add_library(HiPhysics ${PHYSICS_SRC})

set_property(TARGET HiPhysics
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(HiPhysics PUBLIC HiBase)
