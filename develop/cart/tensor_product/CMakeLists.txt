#===============================================================================
# TPO Elliptic CI

add_executable( validate__cart__tpo_elliptic_ci
                validate__cart__tpo_elliptic_ci.f )

target_link_libraries(validate__cart__tpo_elliptic_ci HiCart)

configure_file( validate__cart__tpo_elliptic_ci.prm
                validate__cart__tpo_elliptic_ci.prm COPYONLY )

#===============================================================================
# TPO Spectral

add_executable( validate__cart__tpo_spectral
                validate__cart__tpo_spectral.f )

target_link_libraries(validate__cart__tpo_spectral HiCart)

configure_file( validate__cart__tpo_spectral.prm
                validate__cart__tpo_spectral.prm COPYONLY )

#===============================================================================
# TPO Spectral Iso

add_executable( validate__cart__tpo_spectral_iso
                validate__cart__tpo_spectral_iso.f )

target_link_libraries(validate__cart__tpo_spectral_iso HiCart)

configure_file( validate__cart__tpo_spectral_iso.prm
                validate__cart__tpo_spectral_iso.prm COPYONLY )

#===============================================================================
# TPO Schwarz

add_executable( validate__cart__tpo_schwarz
                validate__cart__tpo_schwarz.f )

target_link_libraries(validate__cart__tpo_schwarz HiCart)

configure_file( validate__cart__tpo_schwarz.prm
                validate__cart__tpo_schwarz.prm COPYONLY )

#===============================================================================
# TPO Schwarz Cubic

add_executable( validate__cart__tpo_schwarz_cubic
                validate__cart__tpo_schwarz_cubic.f )

target_link_libraries(validate__cart__tpo_schwarz_cubic HiCart)

configure_file( validate__cart__tpo_schwarz_cubic.prm
                validate__cart__tpo_schwarz_cubic.prm COPYONLY )
