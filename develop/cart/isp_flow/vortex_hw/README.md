## Vortex_HW

2D Navier-Stokes solution featuring a vortex array proposed in Hesthaven & Warburtin (2008). The example is extended to 3D by imposing periodicity in the z-direction. 

### Temporal convergence

This test explores the convergence of the SDC method. The number of subintervals $$M$$ and correction sweeps $$K$$ are set by specifying the parameters `n_sub` and `n_sweep`  in `isp_flow__sdc_test.prm`. The  test can be run interactively with

    bash -l sdc_m.sh

If enough cores are available it can be run also with 4 processes

    export NP1=2
    export NP2=2
    export EP=2
    bash -l sdc_m.sh

Here,  `NP1` and  `NP2` are the number of partitions in directions 1 and 2, respectively, and `EP` the number of elements in these directions. The convergence results can be visualized with

    python3 convergence_dt.py

The produced graph is stored in `convergence_dt.pdf`.  Here are the reference results for $$M=K=1$$.

![Temporal convergence](convergence_dt.pdf)

On a HPC system you may run larger tests using the provided slurm script, e.g.

    sbatch convergence_dt.slurm

### Spatial convergence

This test studies the $h$-convergence with fixed polynomial degrees $$P_u, P_p$$ for velocity and pressure, respectively, and $P_q$  for integrating the convective term. The polynomial degrees are specified in `isp_flow__sdc_test.prm` by parameters `po_u`, `po_p`  and `po_q`. The time step $$\Delta t$$ is chosen small enough to render the time discretization error negligible.

To limit the runtime it is recommended to reduce the final time $$T$$ by adjusting the corresponding input parameter, e.g. `t_end = 0.1`. The test can be run with

    bash -l convergence_dx.sh

or with slurm

    sbatch convergence_dx.slurm

The results are processed with

    python3 convergence_dx.py

which produces the output file `convergence_dx.pdf`. The reference results shown below were obtained with $$P_u=5$$, $$P_p=4$$  and  $$P_q=9$$.

![Spatial convergence](convergence_dx.pdf)

### References

J.S. Hesthaven & T. Warburton, Nodal Discontinuous Galerkin Methods, Springer 2008