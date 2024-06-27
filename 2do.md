# To Do
[toc]
## HDF5
- [ ] Observe bug report. When fixed remove auxiliary data structures and instructions.
- [ ] I/O of multilevel data

## Adaptive Multilevel Techniques

- [ ] Move `refinement` from `ML_Mesh_3D` into mesh components?!
- [ ] Generalize adapation to accomodate cloning
  - [ ] revise adaptation and refinement marks
  - [ ] generalize existing types and procedures
    - [ ] `ChildDistributionMap_3D` , consider setting `id_child(:,:,:,e)` to identical value
    - [ ] `BuildChildData`, may require separate version
    - [ ] `BuildConnections` etc. 

- [ ] `ML_Mesh_3D` 
  - [ ] Initialization
  - [ ] Generate sequence of clones and/or regular refinements

- [ ] `ML_Adaption_3D`
- [ ] Revise refinement types
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

