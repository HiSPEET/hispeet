## Vortex_TG

Traveling 2D Taylor:Green vortex example proposed in M. Minion & R. Saye, _Higher-order temporal integration for the incompressible Navier-Stokes equations in bounded domains_, J. Comput. Phys. 351: 797-822, 2018. The example is extended to 3D by imposing periodicity in the z-direction, see J. Stiller, _A spectral deferred correction method for incompressible flow with variable viscosity_, [arXiv:2001.11902](https://arxiv.org/abs/2001.11902), 2020. 

### Standalone test

Test case with fixed parameters that is small enough to be run on a single core in interactive mode. It reads the default input file `isp_flow__sdc_test.prm` and the problem-specific file `vortex_tg.prm`. For more details on the parameters see the problem module defined in `isp_flow_problem__vortex_tg.f`.

Run the test with

    mpirun -n 1 ../isp_flow__sdc_test

Upon successful completion it produces the Partitioned VTK Unstructured Data file `vortex_tg.pvtu`, which can be visualized using ParaView. The following example shows the z-component of the computed vorticity at $$t=0.25$$. 

![vorticity](vortex_tg.png)

