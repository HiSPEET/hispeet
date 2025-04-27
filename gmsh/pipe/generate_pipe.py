#> summary:  Mesh for pipe flow
#> author:   Moritz Kreuseler
#> date:     2025/03/21
#> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
#>
#> OH mesh topology:
#>
#>       _ _ _nc_ _ _
#>      /\          /\ 
#>     /  nr       nr \
#>    /    \ _nc_ /    \
#>   /      |    |      \
#>  nc     nc  H nc     nc
#>   \      |_nc_|      /
#>    \    /      \    /
#>     \ nr        nr /
#>      \/_ _ nc _ _\/
#>
#===============================================================================

import gmsh
import os
import sys
import numpy as np
import math
import code

sys.path.append("../utils/")
from gmsh_tools import *

gmsh.initialize()

def createGeometryAndMesh(casename, l, radius, n_c, n_r, n_z, prog):

    #---------------------------------------------------------------------------
    # Initialize
    
    gmsh.clear()
    gmsh.model.add(casename)

    M_array_z, M_z = initAffineTransform(z_transl=-l)

    # to create n elements, gmsh requires the input n+1
    n_c += 1
    n_r += 1
    n_z += 1

    #---------------------------------------------------------------------------
    # Build geometry (from top to bottom ;)

    # contours of first periodic surface: top ..................................

    # cylinder center
    gmsh.model.geo.addPoint( 0, 0, l/2, tag=0)

    # outer quarter points
    gmsh.model.geo.addPoint( 0      , radius   , l/2, tag=1)
    gmsh.model.geo.addPoint(-radius , 0        , l/2, tag=2)
    gmsh.model.geo.addPoint( 0      ,-radius   , l/2, tag=3)
    gmsh.model.geo.addPoint( radius , 0        , l/2, tag=4)

    # inner quarter points on half radius
    gmsh.model.geo.addPoint( 0       , radius/2, l/2, tag=11)
    gmsh.model.geo.addPoint(-radius/2, 0       , l/2, tag=12)
    gmsh.model.geo.addPoint( 0       ,-radius/2, l/2, tag=13)
    gmsh.model.geo.addPoint( radius/2, 0       , l/2, tag=14)

    # outer bound edges: circle 
    gmsh.model.geo.addCircleArc(1, 0, 2, tag=1)
    gmsh.model.geo.addCircleArc(2, 0, 3, tag=2)
    gmsh.model.geo.addCircleArc(3, 0, 4, tag=3)
    gmsh.model.geo.addCircleArc(4, 0, 1, tag=4)

    # inner square h-mesh block
    gmsh.model.geo.addLine(11, 12, tag=11)
    gmsh.model.geo.addLine(12, 13, tag=12)
    gmsh.model.geo.addLine(13, 14, tag=13)
    gmsh.model.geo.addLine(14, 11, tag=14)   

    # quarter separators
    gmsh.model.geo.addLine(1, 11, tag=21)
    gmsh.model.geo.addLine(2, 12, tag=22)
    gmsh.model.geo.addLine(3, 13, tag=23)
    gmsh.model.geo.addLine(4, 14, tag=24)

    gmsh.model.geo.synchronize()

    # contours of bottom and wall ..............................................

    # extrude points
    for i in range(1,5):
        gmsh.model.geo.extrude([(0,i)], 0, 0, -l)
    for i in range(11,15):
        gmsh.model.geo.extrude([(0,i)], 0, 0, -l)

    # extrude outer circle
    for i in range(1,5):
        gmsh.model.geo.extrude([(1,i)], 0, 0, -l)
    surf_wall = [36, 40, 44, 48]

    # extrude inner square h-mesh block
    for i in range(11,15):
        gmsh.model.geo.extrude([(1,i)], 0, 0, -l)
    surf_sqr = [52, 56, 60, 64]

    # extrude quarter separators
    for i in range(21,25):
        gmsh.model.geo.extrude([(1,i)], 0, 0, -l)
    surf_qtr = [68, 72, 76, 80]
    
    # group all 1D entities (curves) for better handling
    curv_cyl    = [ 1,  2,  3,  4, 33, 37, 41, 45] #< nc elements
    curv_sqr    = [11, 12, 13, 14, 49, 53, 57, 61] #< nc elements
    curv_stream = [25, 26, 27, 28, 29, 30, 31, 32] #< nz elements
    curv_radial = [21, 22, 23, 24, 65, 69, 73, 77] #< nr elements

    # periodic surfaces ........................................................
    
    # bottom
    gmsh.model.geo.addCurveLoop([-21,  1, 22,-11], tag=1)
    gmsh.model.geo.addCurveLoop([-22,  2, 23,-12], tag=2)
    gmsh.model.geo.addCurveLoop([-23,  3, 24,-13], tag=3)
    gmsh.model.geo.addCurveLoop([-24,  4, 21,-14], tag=4)
    gmsh.model.geo.addCurveLoop([ 11, 12, 13, 14], tag=5)
    surf_bottom = [1,2,3,4,5]

    # top
    gmsh.model.geo.addCurveLoop([-65, 33, 69,-49], tag=6)
    gmsh.model.geo.addCurveLoop([-69, 37, 73,-53], tag=7)
    gmsh.model.geo.addCurveLoop([-73, 41, 77,-57], tag=8)
    gmsh.model.geo.addCurveLoop([-77, 45, 65,-61], tag=9)
    gmsh.model.geo.addCurveLoop([ 49, 53, 57, 61], tag=10)
    surf_top = [6,7,8,9,10]

    for i in surf_top:
        gmsh.model.geo.add_plane_surface([i], i)
    for i in surf_bottom:
        gmsh.model.geo.add_plane_surface([i], i)

    gmsh.model.geo.synchronize()
    
    #---------------------------------------------------------------------------
    # Build volume

    # create closed surface loop and volume
    gmsh.model.geo.addSurfaceLoop([ 1,  6, 52, 68, 36, 72], tag=1)
    gmsh.model.geo.addSurfaceLoop([ 2,  7, 56, 72, 40, 76], tag=2)
    gmsh.model.geo.addSurfaceLoop([ 3,  8, 60, 76, 44, 80], tag=3)
    gmsh.model.geo.addSurfaceLoop([ 4,  9, 64, 80, 48, 68], tag=4)
    gmsh.model.geo.addSurfaceLoop([ 5, 10, 52, 56, 60, 64], tag=5)
    
    for i in range(1,6):
        gmsh.model.geo.addVolume([i],i)
    volumes = [1,2,3,4,5]

    gmsh.model.geo.synchronize()

    #---------------------------------------------------------------------------
    # Transfine: set number of high order elements on curves

    # 1D
    for i in curv_stream:
        gmsh.model.geo.mesh.setTransfiniteCurve(i, n_z, "Progression", prog)
    for i in curv_cyl:
        gmsh.model.geo.mesh.setTransfiniteCurve(i, n_c, "Progression", prog)
    for i in curv_sqr:
        gmsh.model.geo.mesh.setTransfiniteCurve(i, n_c, "Progression", prog)
    for i in curv_radial:
        gmsh.model.geo.mesh.setTransfiniteCurve(i, n_r, "Progression", prog)
    
    # 2D
    for i in surf_wall:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    for i in surf_bottom:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    for i in surf_top:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    for i in surf_qtr:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    for i in surf_sqr:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)

    # 3D
    for i in volumes:
        gmsh.model.geo.mesh.setTransfiniteVolume(i)

    gmsh.model.geo.synchronize()
    
    #---------------------------------------------------------------------------
    # Recombine: recombine triangles to quadrangles/ pyramides to hexas

    gmsh.model.mesh.createGeometry()

    for i in surf_wall:
        gmsh.model.geo.mesh.setRecombine(2, i)
    for i in surf_bottom:
        gmsh.model.geo.mesh.setRecombine(2, i)
    for i in surf_top:
        gmsh.model.geo.mesh.setRecombine(2, i)
    for i in surf_qtr:
        gmsh.model.geo.mesh.setRecombine(2, i)
    for i in surf_sqr:
        gmsh.model.geo.mesh.setRecombine(2, i)

    for i in volumes:
        gmsh.model.geo.mesh.setRecombine(3, i)
    
    gmsh.model.geo.synchronize()
    
    #---------------------------------------------------------------------------
    # Periodicity
    
    gmsh.model.mesh.setPeriodic(2, surf_top, surf_bottom, M_array_z) 

    #---------------------------------------------------------------------------
    # Physical Groups
    
    gmsh.model.addPhysicalGroup(2, surf_wall,   name=f"wall")
    gmsh.model.addPhysicalGroup(2, surf_top,    name=f"top")
    gmsh.model.addPhysicalGroup(2, surf_bottom, name=f"bottom")
    gmsh.model.addPhysicalGroup(3, volumes,     name=f"fluid")

    #---------------------------------------------------------------------------
    # Mesh Generation
    #---------------------------------------------------------

    #showModel()
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(7)
    
    print(" ++++++++++ SAVING +++++++++++")
    gmsh.write(casename+".msh")
    showModel()

#===============================================================================
# user input

# domain size
length   = 3
radius   = 1

# discretization
n_c  = 2    #< circum, per quarter segment
n_r  = 2    #< radial, must be > 1
n_z  = 3    #< stream, periodic direction, must be > 2
prog = 1    #< element size progression, 1 for equidistant elements

#===============================================================================
# main

# mesh name
casename = "pipe_"+str(n_c)+"-"+str(n_r)+"-"+str(n_z)

if (n_r > 1 and n_z > 2):
    createGeometryAndMesh(casename, length, radius, n_c, n_r, n_z, prog)
else:
    print('+++ Stop: Chose nr > 1 and nz > 2 to create a valid mesh')
