#===============================================================================
# Incompressible Navier-Stokes library

add_library(HiIncompressibleNS 
            ins__operator__3d.f)

set_property(TARGET HiIncompressibleNS
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(HiIncompressibleNS PUBLIC HiBase)
