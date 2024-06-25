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

- Revise refinement types
  - new type _clone_ to allow for $$p$$​- or no refinement
  - treat $$p$$-refinement like $$h$$-refinement using active, frozen and absent elements
  - adapt usage in existing procedures


| type | refinement |  value   | children |  role   |
|:----:|:---------- | --------:|:--------:| ------- |
| `-`  | none       |    `-1`  |    0     |         |
| `p`  | clone      |     `0`  |    1     | active  |
| `h`  | regular    |   `100`  |    8     | active  |
| `h`  | regular    |    `50`  |    8     | closure |
| `h`  | face 1:6   |   `1:6`  |    4     | closure |
| `h`  | edge 1:12  |  `7:18`  |    2     | closure |
| `h`  | vertex 1:8 | `19:26`  |    1     | closure |

- Adaptive multilevel mesh
    - array of meshes
    - auxiliary data: element & interpolation operators?
- Adaptive multilevel spectral element mesh
    - sequence of single-type refinements ($$p$$, $$h$$​, identity)
    - 1:1 correspondence to multilevel mesh
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

