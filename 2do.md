# To Do
[toc]
## HDF5
- [ ] Observe bug report. When fixed remove auxiliary data structures and instructions.
- [ ] I/O of multilevel data

## Adaptive Multilevel Techniques

- [ ] Remove usage of `ParentToChildInterpolation_1D` 
- [ ] Generalize adapation to accomodate cloning
- [ ] `ML_Mesh_3D` 
  - [x] Generate sequence of clones and/or regular refinements
  
- [ ] `ML_Adaption_3D`
- [ ] 

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

