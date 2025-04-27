#> summary:  Mesh for channel flow
#> author:   Benedikt Wex, Moritz Kreuseler
#> date:     2025/03/21
#> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
#===============================================================================

import gmsh
import os
import sys
import numpy as np
import math

sys.path.append("../utils/")
from gmsh_tools import *

gmsh.initialize()

def createGeometryAndMesh(casename, L, N):

    length = L[0] # stream
    width  = L[1] # span
    height = L[2] # wall normal

    nx = N[0] + 1
    ny = N[1] + 1
    nz = N[2] + 1

    lc = 0.05
    gmsh.clear()
    gmsh.model.add(casename)

    MArray_x, M_x = makeAffineTransformationMatrix(x_trans=length)
    MArray_y, M_y = makeAffineTransformationMatrix(y_trans=width)

    print(MArray_y)
    print(MArray_x)

    #---------------------------------------------------------------------------
    # Build first spanwise plane

    gmsh.model.occ.addPoint(0,0,0,           tag=1)
    gmsh.model.occ.addPoint(length,0,0,      tag=2)
    gmsh.model.occ.addPoint(length,0,height, tag=4)
    gmsh.model.occ.addPoint(0,0,height,      tag=5)

    # Note: for the progression towards walls, the direction is essential
    gmsh.model.occ.addLine(1,2, tag=1)
    gmsh.model.occ.addLine(2,4, tag=2)
    gmsh.model.occ.addLine(4,5, tag=5)
    gmsh.model.occ.addLine(5,1, tag=6)

    gmsh.model.occ.synchronize()

    gmsh.model.occ.addWire([1,2, 5, 6], tag=1)
    gmsh.model.occ.addPlaneSurface([1],  tag=1)

    gmsh.model.occ.synchronize()
    
    #---------------------------------------------------------------------------
    # Transform spanwise planes with affine matrix

    F2  = gmsh.model.occ.copy([(2,1)])

    gmsh.model.occ.affineTransform(F2,MArray_y)
    
    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    #---------------------------------------------------------------------------
    # Build in and out planes

    # in
    gmsh.model.occ.addLine( 5, 9, tag=21)
    gmsh.model.occ.addLine( 6, 1, tag=22)
    gmsh.model.occ.addWire([ -6,-21, 10, 22], tag=31)
    gmsh.model.occ.addPlaneSurface([31], tag=31)

    #out
    gmsh.model.occ.addLine( 4,  8, tag=31)
    gmsh.model.occ.addLine( 7,  2, tag=32)
    gmsh.model.occ.addWire([ 2, 31, -8, 32], tag=41)
    gmsh.model.occ.addPlaneSurface([41], tag=41)

    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    #---------------------------------------------------------------------------
    # Build walls

    # top
    gmsh.model.occ.addWire([ 1,-32, -7,  22], tag=51)
    gmsh.model.occ.addPlaneSurface([51], tag=51)
    
    # bottom
    gmsh.model.occ.addWire([ -5, 31,  9, -21], tag=61)
    gmsh.model.occ.addPlaneSurface([61], tag=61)

    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    # group all 1D entities (curves) for better handling
    curv_wall  = [  2,  6,  8, 10]
    curv_span  = [ 21, 22, 31, 32]
    curv_stream = [  1,  5,  7,  9]

    # group all 2D entities (sufaces) for better handling
    surf = [  1,  2,  31, 41, 51, 61]

    #---------------------------------------------------------------------------
    # Build volume

    gmsh.model.occ.addSurfaceLoop(surf, tag=1)
    gmsh.model.occ.addVolume([1], tag=1)
    
    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    
    #---------------------------------------------------------------------------
    # Transfine

    # 1D
    for i in curv_wall:
        gmsh.model.mesh.setTransfiniteCurve(i, nz)#, "Progression", prog)
    for i in curv_stream:
        gmsh.model.mesh.setTransfiniteCurve(i, ny)#, "Progression", 1)
    for i in curv_span:
        gmsh.model.mesh.setTransfiniteCurve(i, nx)#, "Progression", 1)

    # 2D
    for i in surf:
        gmsh.model.mesh.setTransfiniteSurface(i)

    # 3D
    gmsh.model.mesh.setTransfiniteVolume(1)

    gmsh.model.occ.synchronize()

    #---------------------------------------------------------------------------
    # Recombine
    
    gmsh.model.mesh.createGeometry()

    for i in surf:
        gmsh.model.mesh.setRecombine(2, i)

    gmsh.model.mesh.setRecombine(3, 1)
    
    gmsh.model.occ.synchronize()

    #---------------------------------------------------------------------------
    # Periodicity

    gmsh.model.mesh.setPeriodic(2, [ 2], [ 1], MArray_y)
    gmsh.model.mesh.setPeriodic(2, [41], [31], MArray_x)
    
    #---------------------------------------------------------------------------
    # Physical Groups

    gmsh.model.addPhysicalGroup(2, [61], name=f"top")
    gmsh.model.addPhysicalGroup(2, [51], name=f"bottom")
    gmsh.model.addPhysicalGroup(2, [31], name=f"inlet")
    gmsh.model.addPhysicalGroup(2, [41], name=f"outlet")
    gmsh.model.addPhysicalGroup(2, [ 1], name=f"left")
    gmsh.model.addPhysicalGroup(2, [ 2], name=f"right")
    gmsh.model.addPhysicalGroup(3, [ 1], name=f"fluid")

    #---------------------------------------------------------------------------
    # Mesh Generation

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(3)
    gmsh.model.occ.synchronize()

    print(" ++++++++++ SAVING +++++++++++")
    gmsh.write(casename+".msh")
    show()

#===============================================================================
# user input

# domain: 
#     stream  , span  , normal
L = [ 2*np.pi , np.pi , 2    ]  # domain size
N = [ 8       , 8     , 8    ]  # n elements

# mesh name
casename = 'channel_'+str(N[0])+'-'+str(N[1])+'-'+str(N[2])

#===============================================================================
# main

createGeometryAndMesh(casename, L, N)
#gmsh.finalize()
