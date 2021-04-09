add_executable(cg_convdiff_1d
    harmonic_wave_package.f
    cg_convdiff_1d.f
    cg_convdiff_1d__utils.f
    cg_convdiff_1d__imex_euler.f
    cg_convdiff_1d__imex_bdf2.f
    cg_convdiff_1d__imex_bdf3.f
    cg_convdiff_1d__imex_rk.f
    cg_convdiff_1d__imex_euler_sdc.f
    cg_convdiff_1d__imex_tr_sdc.f
    cg_convdiff_1d__imex_rk_sdc.f
    )

target_link_libraries(cg_convdiff_1d HiBase)

#===============================================================================
# Configure files

configure_file(cg_convdiff_1d.prm  cg_convdiff_1d.prm  COPYONLY)
configure_file(convergence_dt.py   convergence_dt.py   COPYONLY)
configure_file(convergence_dt.tmpl convergence_dt.tmpl COPYONLY)
configure_file(convergence_dt.sh   convergence_dt.sh   COPYONLY)

