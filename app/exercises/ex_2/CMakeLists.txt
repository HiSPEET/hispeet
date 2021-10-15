add_executable(sem_basics EXCLUDE_FROM_ALL sem_basics.f)
target_link_libraries(sem_basics HiBase)
configure_file(sem_basics.prm sem_basics.prm COPYONLY)

file( COPY sem_basics.py
      DESTINATION
           ${CMAKE_CURRENT_BINARY_DIR} )

configure_file(README.md README.md COPYONLY)

