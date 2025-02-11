## Variable Viscosity

### Copy from Vortex TG -- to be adapted

Traveling 2D Taylor-Green vortex example proposed in M. Minion & R. Saye, _Higher-order temporal integration for the incompressible Navier-Stokes equations in bounded domains_, J. Comput. Phys. 351: 797-822, 2018. The example is extended to 3D by imposing periodicity in the z-direction, see J. Stiller, _A spectral deferred correction method for incompressible flow with variable viscosity_, [arXiv:2001.11902](https://arxiv.org/abs/2001.11902), 2020. 

### Problem parameters

The following problem parameters are set in the file `vortex_tg.prm`

| name     | description                             |
|----------|-----------------------------------------|
| `stokes` | set `T` for switching to Stokes problem |
| `nu`     | kinematic viscosity $$\nu$$             |
| `vt`     | superimposed translation velocity       |
| `xt`     | initial displacement                    |

All studies are performed in the domain $$\Omega = [-\frac{1}{2},\frac{1}{2}]^2\times [0,l_3]$$ with Dirichlet conditions in directions 1 and 2, and periodic conditions in direction 3. For spatial discretization the domain is decomposed into a regular mesh of cubic elements of size $$\Delta x$$. The extension in direction 3 is automatically adjusted to one element layer, such that $$l_3 = \Delta x$$. Introducing the parameters

- $$N_P$$ -- number of mesh partitions  in directions 1 and 2 and
- $$E_P$$ -- number of elements per partition in both directions

yield  the mesh spacing $$\Delta x = \frac{1}{N_P E_P}$$ . In the run scripts these parameters are denoted by `NP`  and `EP`. 

The solution variables $${\bf u} =(\vec v,p)$$ are approximated using a nodal polynomial basis of the order `PO_U`. Pressure computation and quadrature of convection terms are based on different point sets of order `PO_P` and `PO_Q` that are automatically computed from `PO_U`.

### Temporal convergence study

A temporal convergence study can performed using the command 

```bash
bash -l convergence_dt.sh
```

The following parameters are set in the script:

- `TIME_METHOD` -- time integration method
  + `1`  -- IMEX Euler
  + `2`  -- IMEX BDF2
  + `3`  -- IMEX Runge-Kutta

- `T_END`  -- final time $$T$$ 
- `DT_MAX`  -- maximum time step $$\Delta t_{\max}$$ 
- `ST_MAX`  -- number of time step reductions

The last two parameters define a series of time steps

​	$$\Delta t_s = 2^{-s/2} \Delta t_{\max}$$   for   $$s=$$`0`  to  `ST_MAX` 

which are considered in the study.

Running the script generates the data file `convergence_dt.dat`, which can be processed by

```bash
python convergence_dt.py
```

The results are visualized in a diagram that is stored in `convergence_dt.pdf`.

![Temporal convergence with default parameters](convergence_dt.pdf)

### Spatial convergence study

A spatial convergence study can performed using the command 

```bash
bash -l convergence_dx.sh
```

The following parameters are set in the script:

- `TIME_METHOD` -- time integration method, see above
- `T_END`  -- final time $$T$$ 
- `DT` -- time step $$\Delta t$$ small enough to neglect the temporal error
- `PO_U`  -- polynomial order
- `NP`  -- array containing a series with the number of partitions in directions 1 and 2
- `EP`  -- array containing a series with the number of elements per partition in directions 1 and 2 

The length of the arrays `NP` and `EP`  defines the number of configurations `NC` and must be the same for both. As a consequence, the mesh spacing varies according to $$\Delta x(i) = \frac{1}{N_P(i) E_P(i)}$$ for $$i=0$$ to `NC`$$-1$$ while the polynomial order remains fixed. To adjust the order it can be set as an environment variable before executing the study, e.g.

```bash
export PO_U=7
```

Running the script generates the data file `convergence_dx.dat`, which can be processed by

```bash
python convergence_dx.py
```

The results are visualized in a diagram that is stored in `convergence_dx.pdf`.

![Spatial convergence with default parameters](convergence_dx.pdf)



### Tests on HPC systems

On HPC systems test are typically run in batch mode. For both studies, exemplary Slurm scripts are provided that are configured for the ZIH *Barnard* cluster at TU Dresden. They can be submitted using the `sbatch` command, i.e.

    sbatch convergence_dt-mpi.slurm

or

    sbatch convergence_dx-mpi.slurm

The number of partitions can be adjusted depending on the available ressources. For the temporal convergence study this can be achieved by adjusting the parameters `NP` and `EP` in `convergence_dt-mpi.slurm`. For example, the pairs `NP=1, EP=8` and `NP=4, EP=2` both result in 8 elements in directions 1-2, but 1 and 16 partitions respectively.

For the spatial convergence study the Slurm scrpit switches to HPC configuration arrays `NP` and `EP` that are defined in `convergence_dx.sh`. As `bash` does not support the export of array variables, they need do adjusted directly in the shell script.

Before submitting the jobs take care to load the required *environment modules* and to check the parameters in the Slurm script.
