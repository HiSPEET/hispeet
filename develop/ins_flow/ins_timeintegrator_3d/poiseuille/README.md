# Poiseuille flow

## Periodic case

- running the test case

    - sequential

             $MPIRUN -n 1 ../ins_timeintegrator_3d_test test-poiseuille-periodic
    
    - parallel, using 4 processes

            $MPIRUN -n 4 ../ins_timeintegrator_3d_test test-poiseuille-periodic

## Open boundary case

- works analogously, e.g.

        $MPIRUN -n 4 ../ins_timeintegrator_3d_test test-poiseuille-open

