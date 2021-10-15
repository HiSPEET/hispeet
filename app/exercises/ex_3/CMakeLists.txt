add_executable(helmholtz_1d EXCLUDE_FROM_ALL helmholtz_1d.f)
target_link_libraries(helmholtz_1d HiBase)

configure_file(helmholtz_1d.prm helmholtz_1d.prm COPYONLY)

configure_file(helmholtz_1d.py helmholtz_1d.py COPYONLY)

configure_file(helmholtz_1d_study.py helmholtz_1d_study.py COPYONLY)

configure_file(README.md README.md COPYONLY)
