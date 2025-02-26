# Poiseuille flow

## Periodic case

- running the test case

    - sequential

             mpirun -n 1 ../ins_timeintegrator_3d_test test_periodic
    
    - parallel, using 4 processes

            mpirun -n 4 ../ins_timeintegrator_3d_test test_periodic

## Open boundary case

- works analogously, e.g.

        mpirun -n 4 ../ins_timeintegrator_3d_test test_open

