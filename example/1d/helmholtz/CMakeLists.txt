#===============================================================================
# Executables

# helmholtz_sem_1d-sequential ..................................................

add_executable(helmholtz_sem_1d-sequential helmholtz_sem_1d-sequential.f)

target_link_libraries(helmholtz_sem_1d-sequential hibase)

#===============================================================================
# Configure files

configure_file(helmholtz_sem_1d.prm helmholtz_sem_1d.prm COPYONLY)
