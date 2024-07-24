# To Do
[toc]
## HDF5
- [ ] Observe bug report. When fixed remove auxiliary data structures and instructions.

## Style of OpenMP sections

- [ ] Fit mark up to existing style elements, e.g.
  
     ``````fortran
     !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !$omp master
     code ..
     !$omp end master
     !$omp barrier
     !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ``````

## Revision of Element Operators
- [ ] Settle discrepancy between 1D element operators and mesh operations, e.g.
  - current form for 1D element operator
  
    ``````fortran
    module Coarse_To_Fine_Interpolation__1D
      type CoarseToFineInterpolation_1D
      type CoarseToFineInterpolationOptions_1D
    ``````
  
  - suggested new form indicating operator type
    ``````fortran
    module Coarse_To_Fine_Interpolation_Operator__1D
        type CoarseToFineInterpolationOperator_1D
        type CoarseToFineInterpolationOptions_1D
    ``````
  
## Adaptive Multilevel Techniques

- [ ] `ML_Mesh_3D` 
  - [x] Generalize adapation to accomodate cloning
  - [x] Generate sequence of clones and/or regular refinements
  - [ ] Recursive partitioning
- [x] `ML_Operators_3D`
  - [x] Coarse-to-fine interpolation operator
  - [x] Fine-to-Coarse interpolation and L2 projection operators
- [ ] Tensorproduct operators
  - [ ] Interpolation 8:1
  - [ ] Interpolation/projection 1:8
  - [ ] LIBXSMM versions
  - [ ] Validation

- [ ] `ML_MeshVariable`
  - [ ] Constructor and destructor
  - [ ] HDF5 read and write
  - [ ] Coarse-to-fine interpolation (global, parallel, adaptive)
  - [ ] Fine-to-coarse interpolation/projection
- [ ] Multigrid for elliptic operators
  - [ ] Dependency analysis
  - [ ] Design of data structures and interfaces
  - [ ] Define further steps



## Legacy

Planned re-integration of features available in _HiSPEET_  precursors (`HiSPEET-legacy`)

      > HiSPEET-old/trunk
        - src/differential_operators
        - src/util
        - program/curved
    
      > HiSPEET-old/branches/js-first-flow
        - elliptic solvers
      
      > HiBASE/trunk/src/mesh/
        - mesh/generic_surface_mesh.f
        - bezier
        - triangles

