## Development
### Extension to SVV and IP-H

* `src/cart/elliptic/cart__elliptic_operator_ip__mp_apply_ci.f`
  - extension to hybridizable case
* `src/cart/isp_flow/cart__isp_flow__operators.f`
  - extension to IP-H
- give `IP_ElementOptions1D` instead constructing them?
  e.g. for initialization of `pmg_u` and `pmg_p`
* testing

## Issues
* Accuracy problem when using Intel with `-O3`
* Vortex_HW example converges somewhat slower than with branch `isp_flow-var_diff`, convergence rate with $$P=5$$  only 5 instead of 6

