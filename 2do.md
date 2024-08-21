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
    type StandardOperators_1D
        type StandardOperatorOptions_1D
        ! vs
    type CG_ElementOperators_1D
    type CG_ElementOptions_1D
    ! vs
    type CoarseToFineInterpolation_1D
    type CoarseToFineInterpolationOptions_1D
        ``````
      
  - proposed changes
    ``````fortran
    type StandardElementOperators_1D               ! added 'Element'
    type StandardElementOptions_1D                 ! removed 'Operator'
    type EmbeddedInterpolationOperator_1D          ! added 'Operator'
    type CoarseToFineInterpolationOperator_1D      ! ..
    type FineToCoarseProjectionOperator_1D         ! ..
    ``````
      
  - some alternatives
    ``````fortran
    type CoarseToFineInterpolationOperator_1D
    type CoarseToFineInterpolation_1D
    type CoarseToFineOperator_1D
    type CtFInterpolationOperator_1D
    type HpRefinementOperator_1D                   ! current favorite
    ..
    type HpCoarsensingOperator_1D                  ! ditto
    ``````
- [ ] Element operator data basis 
  - avoid recomputation of element operators
  - optional initialization of standard operators for range of polynomial degrees
  - alternatively or additionally, tabularization of collocation points and weights

## Adaptive Multilevel Techniques

- [ ] `ML_Mesh_3D` 
  - [x] Generalize adapation to accomodate cloning
  - [x] Generate sequence of clones and/or regular refinements
  - [ ] Recursive partitioning
- [x] `ML_Operators_3D`
  - [x] Coarse-to-fine interpolation operator
  - [x] Fine-to-Coarse interpolation and L2 projection operators
- [x] Tensorproduct operators
  - [x] Interpolation 8:1
  - [x] Interpolation/projection 1:8 (TBD: XSMM transition threshold)
  - [x] LIBXSMM versions
- [ ] `ML_MeshVariable`
  - [x] Constructor and destructor
  - [x] HDF5 read and write
  - [x] Parent-to-child interpolation (global, parallel, adaptive)
  - [x] Child-to-parent interpolation/projection
  - [ ] Validation (TBD: HDF5 I/O)
- [ ] Elliptic operators with variable diffusivity
  - [x] Elliptic tensor-product operators for variable diffusivity
  - [x] Extend operator and solvers to variable diffusivity
  - [x] Verification
  - [ ] Serial performance
  - [ ] Scalability with MPI, OpenMP and hybrid
- [ ] Multigrid for elliptic operators
  - [ ] Design of data structures and interfaces
  - [ ] Initialization
  - [ ] FAS multigrid method
  - [ ] FAS-accelerated multilevel Krylov method
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

