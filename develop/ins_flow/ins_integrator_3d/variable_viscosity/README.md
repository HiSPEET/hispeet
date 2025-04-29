## Variable Viscosity

Test case featuring a traveling 3D vortex with variable viscosity. For a detailed description see *J. Stiller, J. Comput. Phys. 423: 109840, 2020* or [arXiv:2001.11902](https://arxiv.org/abs/2001.11902).

### Problem parameters

The following problem parameters are set in the file `variable_viscosity.prm`

| name        | description                               |
| ----------- | ----------------------------------------- |
| `stokes`    | set `T` for switching to Stokes problem   |
| `test_case` | type of the viscosity fluctuation         |
| `nu_0`      | constant baseline viscosity $$\nu_0$$     |
| `nu_1`      | amplitude of variable viscosity $$\nu_1$$ |

Three types of viscosity fluctuations are available:

1. spatially varying viscosity $$\nu(\textbf{x})$$, 
2. spatio-temporally varying viscosity $$\nu(\textbf{x},t)$$,
3. nonlinear viscosity $$\nu(\textbf{u})$$ based on the kinetic energy of the approximate velocity.

The examples use a domain $$\Omega = [-\frac{1}{2},\frac{1}{2}]^3$$ with Dirichlet or periodic conditions . For spatial discretization the domain is decomposed into a regular mesh of cubic elements of size $$\Delta x$$. Choosing the paramaters 

- $$N_{P}$$ -- number of mesh partitions  per direction and
- $$E_{P}$$ -- number of elements per partition in each direction

yields  the mesh spacing $$\Delta x = \frac{1}{N_P E_P}$$ . In the run scripts these parameters are denoted by `NP`  and `EP`. 

The solution variables $${\bf u} =(\vec v,p)$$ are approximated using a nodal polynomial basis of the order `PO_U`. Pressure computation and quadrature of convection terms are based on different point sets of order `PO_P` and `PO_Q` that are automatically computed from `PO_U`. The following test cases are provided:

- `test_1g.prm` -- example using a single-grid solver
- `test_mg-p.prm` -- example using a $p$-multigrid pressure solver
- `test_mg-hp.prm` -- example using an $hp$-multigrid pressure solver

### Temporal convergence study

:warning: *is not yet available*

### Spatial convergence study

 :warning: *is not yet available*

### Temporal stability study

 :warning: *is not yet available*

