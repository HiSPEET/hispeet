#===============================================================================
# TPO Elliptic CI

add_executable( validate__cart__tpo_elliptic_ci
                validate__cart__tpo_elliptic_ci.f )

target_link_libraries(validate__cart__tpo_elliptic_ci HiCart)

configure_file( validate__cart__tpo_elliptic_ci.prm
                validate__cart__tpo_elliptic_ci.prm COPYONLY )

#===============================================================================
# TPO Schwarz

add_executable( validate__cart__tpo_schwarz
                validate__cart__tpo_schwarz.f )

target_link_libraries(validate__cart__tpo_schwarz HiCart)

configure_file( validate__cart__tpo_schwarz.prm
                validate__cart__tpo_schwarz.prm COPYONLY )
