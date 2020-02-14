#===============================================================================
# Executables

# isp_sdc_flow_test ............................................................

add_executable(isp_flow__sdc_test isp_flow__sdc_test.f)

target_link_libraries(isp_flow__sdc_test HiCart)

#===============================================================================
# Configure files

# input file
configure_file(isp_flow__sdc_test.prm isp_flow__sdc_test.prm COPYONLY )

# example case files
set(CASES stokes_gms.prm)
foreach(case ${CASES})
    configure_file(${case} ${case} COPYONLY)
endforeach(case)

# files for batch processing
configure_file( isp_flow__sdc_test.tmpl   isp_flow__sdc_test.tmpl   COPYONLY )
configure_file( isp_flow__sdc_test.sh     isp_flow__sdc_test.sh     COPYONLY )
configure_file( isp_flow__sdc_test.slurm  isp_flow__sdc_test.slurm  COPYONLY )

#===============================================================================
# Test cases

add_subdirectory(stokes_dkm)
add_subdirectory(vortex_hw)
add_subdirectory(vortex_tg)
add_subdirectory(variable_viscosity)

