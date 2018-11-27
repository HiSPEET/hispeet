#===============================================================================
# TPO Grad

add_executable( validate__cart__tpo_grad
                validate__cart__tpo_grad.f )

target_link_libraries(validate__cart__tpo_grad HiCart)

configure_file( validate__cart__tpo_grad.prm
                validate__cart__tpo_grad.prm COPYONLY )

#===============================================================================
# TPO Div

add_executable( validate__cart__tpo_div
                validate__cart__tpo_div.f )

target_link_libraries(validate__cart__tpo_div HiCart)

configure_file( validate__cart__tpo_div.prm
                validate__cart__tpo_div.prm COPYONLY )

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
# TPO Schwarz Iso

add_executable( validate__cart__tpo_schwarz_iso
                validate__cart__tpo_schwarz_iso.f )

target_link_libraries(validate__cart__tpo_schwarz_iso HiCart)

configure_file( validate__cart__tpo_schwarz_iso.prm
                validate__cart__tpo_schwarz_iso.prm COPYONLY )
