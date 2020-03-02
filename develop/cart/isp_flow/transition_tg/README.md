## Transition_TG

3D Taylor Green vortex example proposed ....
Problem is fully periodic.

### Standalone test

Test case with fixed parameters that is small enough to be run on a single core in interactive mode. It reads the default input file `isp_flow__sdc_test.prm` and the problem-specific file `transition_tg.prm`. For more details on the parameters see the problem module defined in `isp_flow_problem__transition_tg.f`.

Run the test with

    mpirun -n 1 ../isp_flow__sdc_test

Upon successful completion it produces the Partitioned VTK Unstructured Data file `transition_tg.pvtu`, which can be visualized using ParaView. 



