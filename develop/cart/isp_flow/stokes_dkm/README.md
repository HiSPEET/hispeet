## Stokes_DKM

This is a 2D Stokes example proposed in M. Deville, L. Kleiser & F. Montigny-Rannou, _Pressure and time treatment for Chebyshev spectral solution of a Stokes problem_, Int. J. Numerical Methods in Fluids, 4(12): 1149-1163, 1984. It is generalized to 3D by imposing periodicity in the z-direction.

The test cases on

  * [temporal convergence](#temporal-convergence)
  * [temporal stability](#temporal-stability)

are adopted from N. Fehn, W.A. Wall & M. Kronbichler, _On the stability of projection methods for the incompressible Navier-Stokes equations based on high-order discontinuous Galerkin discretizations_, J. Comput. Phys. 351: 392-421, 2017.

Befor running the tests, take care that the environment matches the settings at build time.

### Temporal convergence

Check and adjust the settings in the input template `isp_flow__sdc_test.tmpl` and in the run script `convergence_dt.sh`. In interactive mode you can run the test just by typing

    bash -l convergence_dt.sh

Note that this command will start an MPI job with one process. If your system provides at least four computing cores, you can configure the script to use 4 MPI processes as follows

    export NP1=2
    export NP2=2
    export EP=2
    bash -l convergence_dt.sh

Here,  `NP1` and  `NP2` are the number of partitions in directions 1 and 2, respectively, and `EP` the number of elements in these directions. 
On a HPC system with `slurm` you may start a batch job with

    sbatch convergence_dt.slurm

The provided slurm script is configured for `taurus` at ZIH (TU Dresden) and requires that you are a member of the `p_hispeet` project. Otherwise the script needs to be adjusted to meet your HPC ennvironment.

Unless you changed `RANGE_N_SUB` in `convergence_dt.sh` the results can be plotted with

    python3 convergence_dt.py

The produced graph is stored in `convergence_dt.pdf`.  For comparison see the corresponding file in this directory.

![temporal convergence results](convergence_dt.pdf).

### Temporal stability

This test case is intended to check the stability against small time steps. Run the test and plot the results with

    bash -l stability_dt.sh
    python3 stability_dt.py

Alternatively the test can be run with slurm

    sbatch stability_dt.slurm

The produced graph is stored in `stability_dt.pdf`.  For comparison see the corresponding file in this directory.

![temporal stabilty results](stability_dt.pdf).

