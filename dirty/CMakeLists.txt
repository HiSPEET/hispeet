add_executable(cg_helmholtz_1d
    cg_helmholtz_1d.f
    cg_element_operators_1d.f
    cg_spectral_element_utils_1d.f
    cg_condensed_elliptic_solver_1d.f
    )

target_link_libraries(cg_helmholtz_1d HiBase)


add_executable(cg_conv_diff_1d
    cg_conv_diff_1d.f
    cg_conv_diff_1d__utils.f
    cg_element_operators_1d.f
    cg_spectral_element_utils_1d.f
    cg_condensed_elliptic_solver_1d.f
    harmonic_wave_package.f
    cg_convdiff_1d__imex_euler.f
    cg_convdiff_1d__imex_bdf2.f
    )

target_link_libraries(cg_conv_diff_1d HiBase)
