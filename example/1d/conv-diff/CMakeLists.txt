add_executable(cg_convdiff_1d
    spectral_deferred_correction.f
    harmonic_wave_package.f
    cg_convdiff_1d.f
    cg_convdiff_1d__utils.f
    cg_convdiff_1d__imex_euler.f
    cg_convdiff_1d__imex_bdf2.f
    cg_convdiff_1d__imex_bdf3.f
    cg_convdiff_1d__imex_rk.f
    cg_convdiff_1d__imex_euler_sdc.f
    cg_convdiff_1d__vark.f
    )

target_link_libraries(cg_convdiff_1d HiBase)

#===============================================================================
# Configure files

configure_file(cg_convdiff_1d.prm cg_convdiff_1d.prm COPYONLY)

