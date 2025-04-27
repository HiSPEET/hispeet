#> summary:  Mesh for flow past 2D cylinder in a channel with
#> author:   Moritz Kreuseler
#> date:     2025/03/25
#> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
#>
#> Mesh dimensions and topology
#>
#>   - cylinder with "diameter" is placed at x = y = 0
#>   - domain structure and size:
#>
#>  |-----------------------length------------------------|
#>  |----width----|
#>   ______nc_____ ________________nx_____________________  _ _
#>  |\           /|                                       |  |
#>  |  nr     nr  |                                       |  |
#>  |    \   /    |                                       |  |  
#>  nc     O      nc                                      nc |width                 
#>  |    /   \    |                                       |  |
#>  |  nr     nr  |                                       |  |
#>  |/_____nc____\|________________nx_____________________| _|_
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

def createGeometryAndMesh(casename, L, N, prog):

    #---------------------------------------------------------------------------
    # Initialize

    gmsh.clear()
    gmsh.model.add(casename)

    h   = L[0]
    w   = L[1]
    l   = L[2]
    d   = L[3]
    n_z = N[0] + 1
    n_x = N[1] + 1
    n_c = N[2] + 1
    n_r = N[3] + 1

    M_array_z, M_z = initAffineTransform(z_transl=h)

    #---------------------------------------------------------------------------
    # Build geometry

    # contours of first periodic surface bottom ................................

    # cylinder center
    gmsh.model.geo.addPoint(0, 0, 0, tag=0)

    # corner points
    gmsh.model.geo.addPoint( w/2   , w/2, 0, tag=1)
    gmsh.model.geo.addPoint(-w/2   , w/2, 0, tag=2)
    gmsh.model.geo.addPoint(-w/2   ,-w/2, 0, tag=3)
    gmsh.model.geo.addPoint( w/2   ,-w/2, 0, tag=4)
    gmsh.model.geo.addPoint( l-w/2 ,-w/2, 0, tag=5)
    gmsh.model.geo.addPoint( l-w/2 , w/2, 0, tag=6)

    # cylinder points
    a = d / (2 * np.sqrt(2) ) 
    gmsh.model.geo.addPoint( a, a, 0, tag=11)
    gmsh.model.geo.addPoint(-a, a, 0, tag=12)
    gmsh.model.geo.addPoint(-a,-a, 0, tag=13)
    gmsh.model.geo.addPoint( a,-a, 0, tag=14)

    # outer bound edges
    gmsh.model.geo.addLine( 1,  2, tag=1)
    gmsh.model.geo.addLine( 2,  3, tag=2)
    gmsh.model.geo.addLine( 3,  4, tag=3)
    gmsh.model.geo.addLine( 4,  1, tag=4)
    gmsh.model.geo.addLine( 4,  5, tag=5)
    gmsh.model.geo.addLine( 5,  6, tag=6)
    gmsh.model.geo.addLine( 1,  6, tag=7)

    # cylinder edges
    gmsh.model.geo.addCircleArc( 11,  0, 12, tag=11)
    gmsh.model.geo.addCircleArc( 12,  0, 13, tag=12)
    gmsh.model.geo.addCircleArc( 13,  0, 14, tag=13)
    gmsh.model.geo.addCircleArc( 14,  0, 11, tag=14)

    # quarter separators
    gmsh.model.geo.addLine( 1, 11, tag=21)
    gmsh.model.geo.addLine( 2, 12, tag=22)
    gmsh.model.geo.addLine( 3, 13, tag=23)
    gmsh.model.geo.addLine( 4, 14, tag=24)

    gmsh.model.geo.synchronize()

    # contours of top and side walls ...........................................

    # extrude points
    for i in range(1,7):
        gmsh.model.geo.extrude([(0,i)], 0,0,h)
    for i in range(11,15):
        gmsh.model.geo.extrude([(0,i)], 0,0,h)

    # extrude outer bounds
    for i in range(1,8):
        gmsh.model.geo.extrude([(1,i)], 0,0,h)
    surf_side = [38, 46, 54, 62]
    surf_in   = [42]
    surf_int  = [50]    
    surf_out  = [58]

    # extrude cylinder
    for i in range(11,15):
        gmsh.model.geo.extrude([(1,i)], 0,0,h)
    surf_cyl = [66, 70, 74, 78]

    # extrude quarters
    for i in range(21,25):
        gmsh.model.geo.extrude([(1,i)], 0,0,h)
    surf_qtr = [82, 86, 90, 94]

    # group all 1D entities (curves) for better handling
    curv_span   = [25, 26, 27, 28, 29, 30, 31, 32, 33, 34]  #< n_z elements
    curv_radial = [21, 22, 23, 24, 79, 83, 87, 91]          #< n_r elements
    curv_cyl    = [11, 12, 13, 14, 63, 67, 71, 75]          #< n_c elements
    curv_sqr    = [ 1,  2,  3,  4,  6, 35, 39, 43, 47, 55]  #< n_c elements
    curv_stream = [ 5,  7, 51, 59]                          #< n_x elements

    # periodic surfaces ........................................................

    # bottom
    gmsh.model.geo.addCurveLoop([ 11,-22, -1, 21], tag=1)
    gmsh.model.geo.addCurveLoop([ 12,-23, -2, 22], tag=2)
    gmsh.model.geo.addCurveLoop([ 13,-24, -3, 23], tag=3)
    gmsh.model.geo.addCurveLoop([ 14,-21, -4, 24], tag=4)
    gmsh.model.geo.addCurveLoop([- 7,- 4,  5,  6], tag=5)
    surf_bottom = [1,2,3,4,5]

    # top
    gmsh.model.geo.addCurveLoop([ 63, -83, -35, 79], tag= 6)
    gmsh.model.geo.addCurveLoop([ 67, -87, -39, 83], tag= 7)
    gmsh.model.geo.addCurveLoop([ 71, -91, -43, 87], tag= 8)
    gmsh.model.geo.addCurveLoop([ 75, -79, -47, 91], tag= 9)
    gmsh.model.geo.addCurveLoop([-59, -47,  51, 55], tag=10)
    surf_top = [6,7,8,9,10]

    for i in surf_bottom:
        gmsh.model.geo.add_plane_surface([i], i)
    for i in surf_top:
        gmsh.model.geo.add_plane_surface([i], i)
    
    gmsh.model.geo.synchronize()
    
    #---------------------------------------------------------------------------
    # Build volume

    # create closed surface loop and volume
    gmsh.model.geo.addSurfaceLoop([ 1,  6, 66, 82, 38, 86], tag=1)
    gmsh.model.geo.addSurfaceLoop([ 2,  7, 70, 86, 42, 90], tag=2)
    gmsh.model.geo.addSurfaceLoop([ 3,  8, 74, 90, 46, 94], tag=3)
    gmsh.model.geo.addSurfaceLoop([ 4,  9, 78, 94, 50, 82], tag=4)
    gmsh.model.geo.addSurfaceLoop([ 5, 10, 62, 50, 54, 58], tag=5)

    for i in range(1,6):
        gmsh.model.geo.addVolume([i],i)
    volumes = [1,2,3,4,5]
    
    gmsh.model.geo.synchronize()

    #---------------------------------------------------------------------------
    # Transfine: set number of high order elements on curves

    # 1D
    for i in curv_span:
        gmsh.model.geo.mesh.setTransfiniteCurve(i, n_z, "Progression", prog)
    for i in curv_cyl:
        gmsh.model.geo.mesh.setTransfiniteCurve(i, n_c, "Progression", prog)
    for i in curv_sqr:
        gmsh.model.geo.mesh.setTransfiniteCurve(i, n_c, "Progression", prog)
    for i in curv_radial:
        gmsh.model.geo.mesh.setTransfiniteCurve(i, n_r, "Progression", prog)
    for i in curv_stream:
        gmsh.model.geo.mesh.setTransfiniteCurve(i, n_x, "Progression", prog)

    # 2D
    for i in surf_cyl:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    for i in surf_out + surf_in + surf_side + surf_int:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    for i in surf_top:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    for i in surf_bottom:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    for i in surf_qtr:
        gmsh.model.geo.mesh.setTransfiniteSurface(i)

    # 3D
    for i in volumes:
        gmsh.model.geo.mesh.setTransfiniteVolume(i)

    gmsh.model.geo.synchronize()

    #---------------------------------------------------------------------------
    # Recombine: recombine triangles to quadrangles/ pyramides to hexas

    gmsh.model.mesh.createGeometry()

    for i in surf_cyl:
        gmsh.model.geo.mesh.setRecombine(2, i)
    for i in surf_out + surf_in + surf_side + surf_int:
        gmsh.model.geo.mesh.setRecombine(2, i)
    for i in surf_top:
        gmsh.model.geo.mesh.setRecombine(2, i)
    for i in surf_bottom:
        gmsh.model.geo.mesh.setRecombine(2, i)
    for i in surf_qtr:
        gmsh.model.geo.mesh.setRecombine(2, i)

    for i in volumes:
        gmsh.model.geo.mesh.setRecombine(3, i)

    gmsh.model.geo.synchronize()
    
    #---------------------------------------------------------------------------
    # Periodicity

    gmsh.model.mesh.setPeriodic(2, surf_top, surf_bottom, M_array_z)

    #---------------------------------------------------------------------------
    # Physical Groups
    
    gmsh.model.addPhysicalGroup(2, surf_in,              name=f"inlet")
    gmsh.model.addPhysicalGroup(2, surf_out,             name=f"outlet")
    gmsh.model.addPhysicalGroup(2, surf_side + surf_cyl, name=f"wall")
    gmsh.model.addPhysicalGroup(2, surf_top,             name=f"top")
    gmsh.model.addPhysicalGroup(2, surf_bottom,          name=f"bottom")
    gmsh.model.addPhysicalGroup(3, volumes,              name=f"fluid")

    #---------------------------------------------------------------------------
    # Mesh Generation

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(7)

    print(" ++++++++++ SAVING +++++++++++")
    gmsh.write(casename+".msh")
    showModel()

#===============================================================================
# user input

# domain size
height   = 0.15 # span wise direction (z)
width    = 0.4  # orthogonal to flow (y)
length   = 2.2  # length of channel starting at cylinder center (x)
diameter = 0.1  # cylinder diameter

# elements
n_z = 3  # spanwise direction (z) requires >2 for periodicity in HiSPEET
n_x = 14 # streamwise direction flow channel (x)
n_c = 2  # circumferential elements per quarter of cylinder (= width of channel)
n_r = 3  # element layers surounding cylinder

#===============================================================================
# main

# mesh name
casename = "cylinder_2d_"+str(n_c)+"-"+str(n_r)+"-"+str(n_x)+"-"+str(n_z)

L = [height, width, length, diameter]
N = [n_z, n_x, n_c, n_r]
prog = 1
createGeometryAndMesh(casename, L, N, prog)
