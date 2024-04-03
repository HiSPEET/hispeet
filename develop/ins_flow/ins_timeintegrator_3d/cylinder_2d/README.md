# 2D laminar flow past a cylinder

## Background

- benchmark for flow past a cylinder placed slightly asymmetrically in a plane duct
- two test cases
  - steady flow
  - unsteady flow, driven by oscillating inflow
- for detailed description see [1] 

## Status

- both test cases run successfully
- for serious testing, the following features need to be provided
  - evaluation of lift and drag
  - restart option
  - acceleration of solvers
- at the outflow, the implementation relies on open boundary conditions that need further validation, as they extend OBC available in literature

## References

1. M. Schaefer, S. Turek, F. Durst, E. Krause, and R. Rannacher. Benchmark computations of laminar flow around a cylinder. In *Flow Simulation with High-Performance Computers II: DFG Priority Research Programme Results 1993-1995*, pages 547–566. Vieweg+Teubner Verlag, Wiesbaden, 1996.