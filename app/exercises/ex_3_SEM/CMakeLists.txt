add_executable(exercise_SEM_2 EXCLUDE_FROM_ALL exercise_SEM_2.f)
target_link_libraries(exercise_SEM_2 HiBase)

configure_file(input.prm input.prm COPYONLY)

configure_file(plot_exercise_2.py plot_exercise_2.py COPYONLY)

configure_file(plot_exercise_2_study.py plot_exercise_2_study.py COPYONLY)

configure_file(README.md README.md COPYONLY)
