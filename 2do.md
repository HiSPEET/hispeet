## Code: extension to SVV and IP-H

* `src/base/spectral_element/ip_element_operators_1d.`
  - revise documentation of subroutines
     `Get_DiffusionMatrix__w_svv`and
     `Get_EllipticEigensystem__w_svv`

* `src/cart/elliptic/cart__elliptic_operator_ip__mp_apply_ci.f`
  - extension to hybridizable case

* `src/cart/elliptic/cart__schwarz_operator.f`
  - adapt documentation of SchwarzOperator3D to cover the SVV case

* `src/cart/isp_flow/cart__isp_flow__operators.f`
  - extension to IP-H
  - give `IP_ElementOptions1D` instead constructing them,
  e.g. for initialization of `pmg_u` and `pmg_p`
  
* testing

## Funding

1. DFG Anträge: MG für variable Koeffizienten und Bewegte Gitter

   * Ziele
       - lineare Löser für variable Koeffizienten?
       - Anwendung auf Phasenfeldmethoden
       - Ursachen für Probleme bei kurzwellige Koeffizienten finden
       - Lösung durch alternative Techniken: FAS, pCG/Schwarz-Glätter, ... ?
       - bessere Skalierbarkeit durch Kombination p/h-MG
