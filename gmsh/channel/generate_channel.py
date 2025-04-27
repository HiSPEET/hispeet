#> summary:  Mesh for channel flow with refinement toward walls
#> author:   Benedikt Wex, Moritz Kreuseler
#> date:     2025/03/25
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

def createGeometryAndMesh(casename, L, N, prog):

    length = L[0] # stream
    width  = L[1] # span
    height = L[2] # wall normal

    n_x = N[0] + 1
    n_y = N[1] + 1
    n_z = N[2] + 2

    lc = 0.05
    gmsh.clear()
    gmsh.model.add(casename)

    M_array_x, M_x = initAffineTransform(x_transl=length)
    M_array_y, M_y = initAffineTransform(y_transl=width)

    print(M_array_y)
    print(M_array_x)

    #---------------------------------------------------------------------------
    # Build first spanwise plane

    gmsh.model.occ.addPoint(0,0,0,             tag=1)
    gmsh.model.occ.addPoint(length,0,0,        tag=2)
    gmsh.model.occ.addPoint(length,0,height/2, tag=3)
    gmsh.model.occ.addPoint(length,0,height,   tag=4)
    gmsh.model.occ.addPoint(0,0,height,        tag=5)
    gmsh.model.occ.addPoint(0,0,height/2,      tag=6)

    # Note: for the progression towards walls, the direction is essential
    gmsh.model.occ.addLine(1,2, tag=1)
    gmsh.model.occ.addLine(2,3, tag=2)
    gmsh.model.occ.addLine(4,3, tag=3)
    gmsh.model.occ.addLine(3,6, tag=4)
    gmsh.model.occ.addLine(4,5, tag=5)
    gmsh.model.occ.addLine(5,6, tag=6)
    gmsh.model.occ.addLine(1,6, tag=7)

    gmsh.model.occ.synchronize()

    gmsh.model.occ.addWire([1,-2, 4, 7], tag=1)
    gmsh.model.occ.addPlaneSurface([1],  tag=1)

    gmsh.model.occ.addWire([-4, -3, 5, -6], tag=2)
    gmsh.model.occ.addPlaneSurface([2],     tag=2)
    gmsh.model.occ.synchronize()
    
    #---------------------------------------------------------------------------
    # Transform spanwise planes with affine matrix

    F11 = gmsh.model.occ.copy([(2,1)])
    F12 = gmsh.model.occ.copy([(2,2)])

    gmsh.model.occ.affineTransform(F11,M_array_y)
    gmsh.model.occ.affineTransform(F12,M_array_y)

    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    #---------------------------------------------------------------------------
    # Build in and out planes

    # in
    gmsh.model.occ.addLine( 5, 11, tag=21)
    gmsh.model.occ.addLine( 6, 10, tag=22)
    gmsh.model.occ.addLine( 1,  7, tag=23)
    gmsh.model.occ.addWire([ 7, 22, 11,-23], tag=31)
    gmsh.model.occ.addWire([ 6, 21, 12,-22], tag=32)
    gmsh.model.occ.addPlaneSurface([31], tag=31)
    gmsh.model.occ.addPlaneSurface([32], tag=32)

    #out
    gmsh.model.occ.addLine( 4, 12, tag=31)
    gmsh.model.occ.addLine( 3,  9, tag=32)
    gmsh.model.occ.addLine( 2,  8, tag=33)
    gmsh.model.occ.addWire([ 2, 32, -9,-33], tag=41)
    gmsh.model.occ.addWire([ 3, 31,-13,-32], tag=42)
    gmsh.model.occ.addPlaneSurface([41], tag=41)
    gmsh.model.occ.addPlaneSurface([42], tag=42)

    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    #---------------------------------------------------------------------------
    # Build walls

    # top
    gmsh.model.occ.addWire([ 5, 31, 14,-21], tag=51)
    gmsh.model.occ.addPlaneSurface([51], tag=51)

    # bottom
    gmsh.model.occ.addWire([ 1, 33, -8,-23], tag=61)
    gmsh.model.occ.addPlaneSurface([61], tag=61)

    # mid
    gmsh.model.occ.addWire([ -4, 32, 10,-22], tag=71)
    gmsh.model.occ.addPlaneSurface([71], tag=71)

    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    # group all 1D entities (curves) for better handling
    curv_wall1  = [  2,  7,  9, 11]
    curv_wall2  = [  3,  6, 12, 13]
    curv_span   = [ 21, 22, 23, 31, 32, 33]
    curv_stream = [  1,  4,  5,  8, 10, 14]

    # group all 2D entities (sufaces) for better handling
    surf = [  1,  2,  3,  4, 31, 32, 41, 42, 51, 61, 71]

    #---------------------------------------------------------------------------
    # Build volume

    gmsh.model.occ.addSurfaceLoop([ 1, 41,  3, 31, 61, 71], tag=1)
    gmsh.model.occ.addSurfaceLoop([ 2, 42,  4, 32, 51, 71], tag=2)

    for i in range(1,3):
        gmsh.model.occ.addVolume([i], tag=i)
    volumes = [1,2]
    
    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    
    #---------------------------------------------------------------------------
    # Transfine

    # 1D
    for i in curv_wall1:
        gmsh.model.mesh.setTransfiniteCurve(i, int(n_z/2), "Progression", prog)
    for i in curv_wall2:
        gmsh.model.mesh.setTransfiniteCurve(i, int(n_z/2), "Progression", prog)
    for i in curv_stream:
        gmsh.model.mesh.setTransfiniteCurve(i, n_y       , "Progression", 1)
    for i in curv_span:
        gmsh.model.mesh.setTransfiniteCurve(i, n_x       , "Progression", 1)

    # 2D
    for i in surf:
        gmsh.model.mesh.setTransfiniteSurface(i)

    # 3D
    for i in volumes:
        gmsh.model.mesh.setTransfiniteVolume(i)

    gmsh.model.occ.synchronize()

    #---------------------------------------------------------------------------
    # Recombine
    
    gmsh.model.mesh.createGeometry()

    for i in surf:
        gmsh.model.mesh.setRecombine(2, i)

    for i in volumes:
        gmsh.model.mesh.setRecombine(3, i)
    
    gmsh.model.occ.synchronize()

    #---------------------------------------------------------------------------
    # Periodicity

    gmsh.model.mesh.setPeriodic(2, [ 3,  4], [ 1,  2], M_array_y)
    gmsh.model.mesh.setPeriodic(2, [41, 42], [31, 32], M_array_x)
    
    #---------------------------------------------------------------------------
    # Physical Groups

    gmsh.model.addPhysicalGroup(2, [51],     name=f"top")
    gmsh.model.addPhysicalGroup(2, [61],     name=f"bottom")
    gmsh.model.addPhysicalGroup(2, [31, 32], name=f"inlet")
    gmsh.model.addPhysicalGroup(2, [41, 42], name=f"outlet")
    gmsh.model.addPhysicalGroup(2, [ 1,  2], name=f"left")
    gmsh.model.addPhysicalGroup(2, [ 3,  4], name=f"right")
    gmsh.model.addPhysicalGroup(3, [ 1,  2], name=f"fluid")

    #---------------------------------------------------------------------------
    # Mesh Generation

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(3)
    gmsh.model.occ.synchronize()

    print(" ++++++++++ SAVING +++++++++++")
    gmsh.write(casename+".msh")
    showModel()

#===============================================================================
# user input

# domain: 
#     stream  , span  , normal
L = [ 2*np.pi , np.pi , 2    ]  # domain size
N = [ 10      , 8     , 12   ]  # n elements

# element width progression towards wall
prog = 1.3

#===============================================================================
# main

# mesh name
casename = 'channel_'+str(N[0])+'-'+str(N[1])+'-'+str(N[2])+'_'+str(prog)

createGeometryAndMesh(casename, L, N, prog)
