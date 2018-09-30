#===============================================================================
# isp_flow

file(GLOB ISP_FLOW isp_flow/*.f)

#===============================================================================
# physics library

list(APPEND PHYSICS_SRC ${ISP_FLOW})

add_library(HiPhysics ${PHYSICS_SRC})

set_property(TARGET HiPhysics
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(HiPhysics PUBLIC HiBase)