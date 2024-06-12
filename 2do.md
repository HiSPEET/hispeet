# To Do
[toc]
## HDF5 Read & Write
### Concept
- HDF5 group structure
    - single level mesh
    - multilevel mesh
    - goal: design SL mesh group to fit in ML group
- Files
    - read/write one file per process
    - detect and handle case where number of files exeeds number of processes
- Writing  (allready implemented in `mesh%WriteHDF5`)
- Reading  (allready implemented in `mesh%ReadHDF5`)

## Flow Solver Restart
### Concept
- Writing
    - HDF5 file for mesh
    - flow data: __TBD__
- Reading
    - read HDF5 mesh file
    - build SEM mesh as usual
    - flow data: __TBD__

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

