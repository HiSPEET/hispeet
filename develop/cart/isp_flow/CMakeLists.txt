#===============================================================================
# Executables

# elliptic_test ................................................................

add_executable(isp_flow__sdc_test isp_flow__sdc_test.f)

target_link_libraries(isp_flow__sdc_test HiCart)

#===============================================================================
# Configure files

# input file
configure_file(isp_flow__sdc_test.prm isp_flow__sdc_test.prm COPYONLY )

# example case files
set(CASES stokes_dkm.prm stokes_gms.prm vortex_hw.prm vortex_tg.prm)
foreach(case ${CASES})
    configure_file(${case} ${case} COPYONLY)
endforeach(case)

# files for batch processing
configure_file( isp_flow__sdc_test.tmpl   isp_flow__sdc_test.tmpl   COPYONLY )
configure_file( isp_flow__sdc_test.sh     isp_flow__sdc_test.sh     COPYONLY )
configure_file( isp_flow__sdc_test.slurm  isp_flow__sdc_test.slurm  COPYONLY )
