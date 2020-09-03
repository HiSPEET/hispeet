## Development
### Extension to SVV and IP-H

* `src/cart/elliptic/cart__elliptic_operator_ip__mp_apply_ci.f`
    - extension to hybridizable case
* `src/cart/isp_flow/cart__isp_flow__operators.f`
    - extension to IP-H
    - give `IP_ElementOptions1D` instead constructing them?
      e.g. for initialization of `pmg_u` and `pmg_p`
* testing

### OpenMP

* Re-integration of OpenMP_Binding module
* Make more routines working with OpenMP
    - `src/cart/elliptic/cart__dg_elliptic_ciu_bc.f`
    - `src/cart/elliptic/cart__elliptic_operator_ip__mp_bctorhs.f`

## Issues

* Accuracy problem when using Intel with `-O3`
* `MPI_REAL_RHP` deactivated since quad precision not yet supported with Intel MPI
* Vortex_HW example converges somewhat slower than with branch `isp_flow-var_diff`, convergence rate with $$P=5$$  only 5 instead of 6

## Delayed upgrades

* Change `add_definitions` to `add_compile_definitions` in CMake files
  - CMakeLists.txt
  - external/cmake/libxsmm.cmake
  etc.

