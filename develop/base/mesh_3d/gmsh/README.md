# Notes on mesh examples

## pipe

- geometry and block structure OK
- try different mesh densities
- try to remove periodicity in axial direction

## cylinder

- upstream domain extension possibly to small
  - compare with test cases defined [1]
  - add block for undisturbed inflow
- try to add periodicity in vertical direction

## naca0012, naca0012_aoa

- block structure needs to be improved to avoid distortions 

## References

1. M. Schaefer, S. Turek, F. Durst, E. Krause, and R. Rannacher. Benchmark computations of laminar flow around a cylinder. In *Flow Simulation with High-Performance Computers II: DFG Priority Research Programme Results 1993-1995*, pages 547–566. Vieweg+Teubner Verlag, Wiesbaden, 1996.