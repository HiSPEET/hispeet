## Development
### Re-Integration

      ~/Numerics/HiSPEET-legacy/HiSPEET-old/trunk
        src/mesh
        - gekrümmtes 3D Gitter (Startvorlage, unvollständig portiert)
        - Integritätsprüfung
        src/differential_operators
        - metrische Koeffizienten, Grad, Div
        src/util
        - Generierung generischer zylindrischer Gitter
        - Import generischer und Exodus-Gitter
        program/curved
        - Moving-Mesh
        
      ~/Numerics/HiSPEET-legacy/HiSPEET-old/branches/js-first-flow
        - elliptische Löser?!
      
      ~/Numerics/HiSPEET-legacy/HiBASE/trunk/src/mesh/
        - generic_volume_mesh.F
        - generic_volume_mesh__import_exodus.f

### Extension to SVV and IP-H

* `src/cart/elliptic/cart__elliptic_operator_ip__mp_apply_ci.f`
    - extension to hybridizable case
* `src/cart/isp_flow/cart__isp_flow__operators.f`
    - extension to IP-H
    - give `IP_ElementOptions1D` instead constructing them?
      e.g. for initialization of `pmg_u` and `pmg_p`
* testing

## Issues

* `MPI_REAL_RHP` deactivated since quad precision not yet supported with Intel MPI
* Vortex_HW example converges somewhat slower than with branch `isp_flow-var_diff`, convergence rate with $$P=5$$  only 5 instead of 6

