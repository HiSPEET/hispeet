MeshElementVertex_3D
  id          ←  IdentifyVertices  ←  element%{edge,face}%{n_neighbor,i_neighbor}
✓ n_neighbor  ←  *** ab initio ***
✓ i_neighbor  ←  *** ab initio ***
  rank        ←  IdentifyRanks     ←  IdentifyVertices
  val         ←  IdentifyRanks

MeshElementEdge_3D
  id          ←  IdentifyEdges     ←  element%vertex%id, element%face%boundary
  orientation ←  IdentifyEdges
✓ n_neighbor  ←  *** ab initio ***
✓ i_neighbor  ←  *** ab initio ***
  rank        ←  IdentifyRanks     ←  IdentifyEdges
  val         ←  IdentifyRanks

MeshElementFace_3D
  id          ←  BuildFaces        ←  element%{vertex,edge}%id
✓ boundary    ←  *** ab initio ***
  normal      ←  BuildFaces
  rotation    ←  BuildFaces
✓ n_neighbor  ←  *** ab initio ***
✓ i_neighbor  ←  *** ab initio ***
  rank        ←  IdentifyRanks     ←  BuildFaces
  val         ←  IdentifyRanks

MeshElement_3D
✓ id          ←  *** ab initio ***
  vertex      →  see MeshElementVertex_3D
  edge        →  see MeshElementEdge_3D
  face        →  see MeshElementFace_3D
✓ neighbor    ←  *** ab initio ***
* adaptation  ←  *** retained ***
✓ geometry    ←  *** ab initio ***











