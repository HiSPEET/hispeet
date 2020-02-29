## Spatial convergence with nonlinear viscosity

This test repeats the $$h$$-convergence study for incompressible flow with a solution-dependent viscosity presented in J. Stiller, _A spectral deferred correction method for incompressible flow with variable viscosity_, [arXiv:2001.11902](https://arxiv.org/abs/2001.11902), 2020. It is rather expensive and, hence, taylored to be run in batch mode on an HPC system.

Before starting the test, inspect the parameters set in the slurm scripts. As a sanity check you may start one job manually, e.g.

```bash
sbatch p06_n04.slurm
```

Inspect the output file and if all works fine terminate the job with `scancel`.  Now you are ready to start the full test

```bash
bash convergence_dx.sh
```

This script will start a series of jobs. Once these are completed, you may evaluate the results as follows

```bash
# generate data from log files
bash log2dat.sh
# generate error plot
python3 convergence_dx.py
```

The Python script generates an error plot and saves it to `convergence_dx.pdf`. If the test succeeded, the plot should look like this:

![Spatial convergence with nonlinear viscosity](convergence_dx.pdf)

