#===============================================================================
# Executables

add_executable( validate__cart__tpo_elliptic_ci
                validate__cart__tpo_elliptic_ci.f )

target_link_libraries(validate__cart__tpo_elliptic_ci HiCart)

#===============================================================================
# Configure files

configure_file( validate__cart__tpo_elliptic_ci.prm
                validate__cart__tpo_elliptic_ci.prm COPYONLY )
