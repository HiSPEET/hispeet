# Generating a GMSH mesh for HiSPEET
 2025/03/31 M. Kreuseler 

[toc]

## 1  Preliminaries

### 1.1  Get GMSH for your python

using pip:

`pip install gmsh`

using conda:

`conda install conda-forge::python-gmsh`

### 1.2  This are the three steps towards a valid HiSPEET mesh

1. Chose the right CAD facility in GMSH. There are two: **"geo"** and **"occ"**
2. Create a geometry by
   - building it using the GMSH CAD kernel **"geo"**
   - building it using the GMSH CAD kernel **"occ"**
   - loading a CAD file and modifying it in GMSH
   - loading a pre-build mesh 
3. Mesh the geometry

## 2  Chose the CAD kernel 
### 2.1  The "occ" kernel
good, because:
  - copying and affine transformation of surfaces is possible
  - import of CAD files or pre-build meshes (e.g. CGNS)
  - periodicity can be defined via an affine transformation matrix (required by HiSPEET)

not so good, because:
 - sometimes produces corrupted meshes for self-build (in GMSH) geometries

### 2.2  The "geo" kernel
good, because:
- securely produces a valid mesh
- preferable for simple self-build geometries
- periodicity can be defined via an affine transformation matrix (required by HiSPEET)

not so good, because:
- weak mesh and CAD import
- surfaces can not be copied and transformed using an affine transformation 

## 3  Build a geometry

### 3.1  Using the GMSH CAD kernel  "geo"

**POINTS: start building the base by adding points  (dim=0 entities)**
`gmsh.model.geo.addPoint( x , y , z , tag=pointA)`
> NOTE
> - `tag` is any integer
> - each entity (e.g. point, line, surface, ...) has a tag that is unique among all entities of the same dimension 

**LINES: connect points to lines (dim=1 entities)**
 `gmsh.model.geo.addLine(pointA, pointB, tag=lineA)`
  - the orientation of line (from `pointA` to `pointB`) is essential for the direction of element refinement using the `progression` function (e.g. for refinement towards walls)
  - the direction of the progression is defined by the orientation of the line: for growing elements towards the left in the example below (fine elements near `pointA` and `pointD`), `lineA` and `lineC` must be oriented from `pointA` to `pointB` and `pointD` to `pointC`, respectively, as indicated by the arrows:

```
			    lineA
		 pointB o-----<-----o pointA
			|           |
		  lineB ^           ^ lineD
			|           |
		 pointC o-----<-----o pointD	 	
		            lineC
				
```

**SURFACES: Build a surface (dim=2 entities) from a loop of dim=1 entities**

 `gmsh.model.geo.addCurveLoop([lineA, -lineB, -lineC, lineD], tag=1DloopA)`

 `gmsh.model.geo.add_plane_surface([1DloopA], tag=surfaceA)` 
 - to create a closed loop, take care on the orientation! 
 - reverse a line with a minus `-` if needed, e.g for the example above to create a closed loop `lineB` and `lineC`must be reversed

**SURFACES: Build a surface (dim=2 entity) by extruding dim=1 entities**

 `gmsh.model.geo.extrude([(0,pointA),(0,pointB),...], dx, dy, dz)`

 `gmsh.model.geo.extrude([(1,lineA) ,(1,lineB) ,...], dx, dy, dz)` 
 - preferable for curved surfaces (e.g. cylinders)
 - always extrude both, relevant points and lines

> NOTE
> A combination `(integer, tag)` (e.g. `(0, pointA)` or `(1, lineB)`)  generally denotes a pair of a "dimension" and a "tag" identifying one specific entity of the model, e.g.:
> `(0, pointA)`  is the **point** (dim=0) with **tag** `pointA`

**VOLUMES: build a volume (dim=3 entity) from a loop of dim=2 entities**

`gmsh.model.geo.addSurfaceLoop([surfaceA, ...], tag=2DloopA)`

`gmsh.model.geo.addVolume([2DloopA], volumeA)`

### 3.2  Using the GMSH CAD kernel : "occ" kernel

 - mostly similar to **"geo"**
 - replace `.geo.` by `.occ.` in all example command lines above
 - additionally, a (periodic) surface can be created by copying an existing surface and applying an affine transformation

 `surfaceB = gmsh.model.occ.copy([(2,surfaceA)])`

 `gmsh.model.occ.affineTransform(surfaceB , T_AB)`
 >NOTE 
 > -  `T_AB` is the 1x16 array representation of a 4x4 affine transformation matrix, a detailed explanation of `T_AB` follows below
 > - the use of a copied, transformed surface sometimes produces a corrupted mesh! 
 > - extruding and affine transformation can not be combined in the same direction ! 

### 3.3  Load and modify a geometry using GMSH ("occ" kernel)
 - this is an advanced topic beyond the scope of this minimum introduction

### 3.4  Load a pre-build mesh ("occ" kernel)

- load a cgns file:

`gmsh.merge("path/to/myFile.cgns")`

`gmsh.model.occ.synchronize()`
> NOTE
> - requirements for the cgns are highly case specific
> - the creation of the pre-build mesh demands careful consideration

## 4  Meshing 

### 4.1  Transfinite (for in-build geometry "geo"/"occ" only) 

- to create a high order hexa mesh, all entities that are supposed to be meshed (dim=1, dim=2, dim=3) must be set as **transfinite**:
	-  set dim=1 entities transfinite: 
	
		`gmsh.model.<geo/occ>.mesh.setTransfiniteCurve(lineA, n, "Progression", r)`
		- `r`: growth parameter of element width (progression) towards one end of the line 
		> NOTE
		> - GMSH uses a geometric progression
		> - the width $\Delta x_i$ of the $i=[1 \ldots n-1]$-th element is $r^{i-1}\Delta x_1$,
		 starting with a element width of $\Delta x_1$ at the first element on the line
		
		- `n`-1: number of elements along this line 
		- the number `n` on opposing lines must match:
```
   			   lineA: n=5
		 pointB o--x---x---x--o pointA
			|	      |
		   	x	      x
  	     lineC: n=4 |             | lineD: n=4
		   	x	      x
			|             |
		 pointC	o--x---x---x--o pointD
			  lineC: n=5
```

  - set dim=2 entities transfinite: 

`gmsh.model.<geo/occ>.mesh.setTransfiniteSurface(surfaceA)`
  - set dim=3 entities transfinite:

`gmsh.model.<geo/occ>.mesh.setTransfiniteVolume(volumeA) `

### 4.2  Recombine (for in-build geometry "geo"/"occ" only) 

- triangles and tetrahedrons on all entities (dim=2 and dim=3) that are set to **transfinite** must be **recombined** to quadrangles and hexas: 

`gmsh.model.<geo/occ>.mesh.setRecombine(2, surfaceA)` 

`gmsh.model.<geo/occ>.mesh.setRecombine(3, volumeA)`

### 4.3  Periodicity
- periodicity is defined, based on the tags an spatial relation of periodic surfaces (dim=2):
  `gmsh.model.mesh.setPeriodic(2, surfaceB, surfaceA, T_AB)`
	- the **affine transformation** between them must be known and given in the 1x16 array representation 
	
		`T_AB` = $$[\mathbf{T}_{11},\mathbf{T}_{12},...,\mathbf{T}_{14},\mathbf{T}_{21},...,\mathbf{T}_{24},\mathbf{T}_{44}]$$ 

      of the 4x4 transformation matrix:
$$
\mathbf{T} = \mathbf{T}_{\mathrm{scale}} \cdot \mathbf{T}_{\mathrm{translate}} \cdot \mathbf{T}_{\mathrm{rot},x} \cdot \mathbf{T}_{\mathrm{rot},y} \cdot \mathbf{T}_{\mathrm{rot},z}
$$
with 
$$
\mathbf{T}_{\mathrm{scale}}  = 
\begin{pmatrix}
sx & 0 & 0 & 0 \\
sy & 0 & 0 & 0 \\
sz & 0 & 0 & 0 \\
0  & 0 & 0 & 1
\end{pmatrix}
$$
$$
\mathbf{T}_{\mathrm{translate}}  = 
\begin{pmatrix}
1 & 0 & 0 & dx \\
0 & 1 & 0 & dy \\
0 & 0 & 1 & dz \\
0  & 0 & 0 & 1
\end{pmatrix}
$$
$$
\mathbf{T}_{\mathrm{rot},x}  = 
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & \cos(\alpha) & \sin(\alpha) & 0 \\
0 & -\sin(\alpha) & \cos(\alpha) & 0 \\
0  & 0 & 0 & 1
\end{pmatrix}
$$
$$
\ldots
$$

>NOTE! In the **periodic direction**, a **minimum of three elements** is required in HiSPEET. An element cannot have the same neighbor on both sides.

### 3.4  Physical groups
- all boundary and periodic surfaces must be defined in physical groups 
- surfaces that are not member of a physical group will be treated as inner surfaces, no boundary conditions ca be imposed there 
- the order of their definition is the order that HiSPEET will read them to the `bc_v` structure, **be careful with this!**

`gmsh.model.addPhysicalGroup(2, [surfaceA], name=f"OUTLET")`

`gmsh.model.addPhysicalGroup(2, [surfaceB], name=f"INLET")`

`gmsh.model.addPhysicalGroup(2, [surfaceC,surfaceD],name=f"BOUND1")`

... 

`gmsh.model.addPhysicalGroup(3, [volumeA,volumeB,...],name=f"FLUID")` 

### 3.5  Mesh generation
- mesh dim=3 entities:

`gmsh.model.mesh.generate(3)`
- set order to 3:

`gmsh.model.mesh.setOrder(3)`
- save the mesh to a `.msh` file:

`gmsh.write("myMesh.msh")` 

## 5  General advice 
- to view the current mesh/geometry in the GMSH GUI define and call this function: 

```
def show():
    def checkForEvent():
        action = gmsh.onelab.getString("ONELAB/Action")
        if len(action) and action[0] == "check":
            gmsh.onelab.setString("ONELAB/Action", [""])
            createGeometryAndMesh()
            gmsh.graphics.draw()
        return True

    if "-nopopup" not in sys.argv:
        gmsh.fltk.initialize()
        while gmsh.fltk.isAvailable() and checkForEvent():
            gmsh.fltk.wait()
    gmsh.finalize()
```
 - if some features are not shown or found correctly, it is benefitial to update the current state by calling:

 `gmsh.model.<geo/occ>.synchronize()`

