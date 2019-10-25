add_executable(exercise_SEM_1 EXCLUDE_FROM_ALL exercise_SEM_1.f)
target_link_libraries(exercise_SEM_1 HiBase)

configure_file(input.prm input.prm COPYONLY)

configure_file(plot_exercise_1.py plot_exercise_1.py COPYONLY)

configure_file(README.md README.md COPYONLY)
