configure_file(input_sem_1d input_sem_1d COPYONLY)

#configure_file(output_1 output_1 COPYONLY)

#configure_file(output_2 output_2 COPYONLY)

add_executable(sem_1d EXCLUDE_FROM_ALL sem_1d.f)
target_link_libraries(sem_1d HiBase)

add_executable(ping_pong EXCLUDE_FROM_ALL ping_pong.f)
target_link_libraries(ping_pong HiBase)

#configure_file(README.md README.md COPYONLY)

