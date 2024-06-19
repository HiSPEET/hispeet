# To Do
[toc]
## HDF5
- Observe bug report. When fixed
    - test
    - remove redundant statements (deallocate) in `mesh%ReadHDF5` 
- Understand HDF5 `groups`
- I/O of multilevel data

## Flow Solver Restart
- Activate and test when HDF5 is fixed

## Adaptive Multilevel Techniques

- Adaptive multilevel mesh
    - datastructures and methods exist
    - final design: *TBD*
- Adaptive multilevel spectral element mesh: __TBD__
    - based on multilevel mesh
    - additional $$p$$-levels
- Interweaving spatial and temporal adaptivity: __TBD__
    - working hypothesis: time mesh can be uniquely mapped to ML spectral element mesh
    - each temporal level maps to one spatial level
    - ML spectral element mesh can be used alone, e.g. for multigrid correction schemes

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

