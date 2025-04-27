import gmsh
import os
import sys
import numpy as np
import math

def initAffineTransform(x_transl=0, y_transl=0, z_transl=0, x_rot=0, y_rot=0, z_rot=0):
    # --------------------------------------------------------
    # Returns the Affine Transform Matrix in Array and Matrix Form for
    # any Translation x,y,z
    # --------------------------------------------------------

    M_scale = np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])

    M_transl = np.array([[1, 0, 0, x_transl],
                         [0, 1, 0, y_transl],
                         [0, 0, 1, z_transl],
                         [0, 0, 0, 1]])

    phi =  x_rot
    R_x  = np.array([[1, 0, 0, 0],
                     [0,  math.cos(phi), math.sin(phi), 0],
                     [0, -math.sin(phi), math.cos(phi), 0],
                     [0, 0, 0, 1]])
 
    phi =  y_rot
    R_y = np.array([[ math.cos(phi), 0, math.sin(phi), 0],
                    [0, 1, 0, 0],
                    [-math.sin(phi), 0, math.cos(phi), 0],
                    [0, 0, 0, 1]])

    phi =  z_rot
    R_z  = np.array([[ math.cos(phi), math.sin(phi), 0, 0],
                     [-math.sin(phi), math.cos(phi), 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])


    M_unity = np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])


    M = M_unity @ M_scale @ R_x @ R_y @ R_z @ M_transl

    M_array = []
    for i in range(M.shape[1]):
        for j in range(M.shape[0]):
            M_array.append(M[i, j])

    return (M_array, M)

def getPointCoord(tag):
    # -------------------------------------------------------------------------
    # Returns coordinates of the point with "tag"
    # -------------------------------------------------------------------------

    gmsh.model.occ.synchronize()

    bounding_tags = gmsh.model.occ.getBoundingBox(0, tag)
    tag_x = bounding_tags[0]+(bounding_tags[3]-bounding_tags[0] )/2
    tag_y = bounding_tags[1]+(bounding_tags[4]-bounding_tags[1] )/2
    tag_z = bounding_tags[2]+(bounding_tags[5]-bounding_tags[2] )/2

    coords = [tag_x, tag_y, tag_z]

    return coords

def getLineLength(dimtag):
    # -------------------------------------------------------------------------
    # Returns length of line of Line defined by dimtag=[dim,tag]
    # -------------------------------------------------------------------------

    gmsh.model.occ.synchronize()

    points = gmsh.model.getBoundary([dimtag])
    values = []

    for point in points:
        value = gmsh.model.getValue(point[0], point[1], [])
        values.append(value)

    start_coord = np.array(values[ 0])
    end_coord   = np.array(values[-1])

    l = np.linalg.norm(end_coord-start_coord)
        
    return(l)

def getLineVector(dimtag):
    # -------------------------------------------------------------------------
    # Returns the vector defining the line dimtag=[dim,tag]
    # -------------------------------------------------------------------------

    gmsh.model.occ.synchronize()

    points = gmsh.model.getBoundary([dimtag])
    values = []

    for point in points:
        value = gmsh.model.getValue(point[0], point[1], [])
        values.append(value)

    start_coord = np.array(values[0])
    end_coord   = np.array(values[-1])

    v = end_coord-start_coord
    
    return v

def importSTL(geompath, scale=1):
    # ------------------------------------------------------------------------
    # Import stl
    # ------------------------------------------------------------------------

    gmsh.model.occ.synchronize()

    # path to .igs geometry
    print('Load geometry from: ' + geompath)

    # load and resize geometry
    gmsh.option.setString('Geometry.OCCTargetUnit', '')
    gmsh.option.setNumber('Geometry.OCCScaling', scale)
 
    gmsh.model.occ.importShapes(os.path.join(geompath), False)

    gmsh.model.occ.synchronize()

def showModel():
    # -------------------------------------------------------------------------
    # Function to show intermediate state at any moment: show()
    # -------------------------------------------------------------------------

    gmsh.model.occ.synchronize()

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


def rotateMesh(x_rot,y_rot,z_rot):
    # ------------------------------------------------------------------------
    # Rotate entire mesh
    # -------------------------------------------------------------------------

    gmsh.model.occ.synchronize()

    # get cartesian coordinates of all nodes    
    node_tags, node_coords,_ = gmsh.model.mesh.getNodes()
    node_coords              = np.array(node_coords).reshape(-1,3)
    
    # build rotation matrices
    R_x  = np.array([[1, 0, 0],
                     [0, math.cos(x_rot), -math.sin(x_rot)],
                     [0, math.sin(x_rot),  math.cos(x_rot)]])

    R_y = np.array([[ math.cos(y_rot), 0, math.sin(y_rot)],
                    [0, 1, 0],
                    [-math.sin(y_rot), 0, math.cos(y_rot)]])

    R_z = np.array([[math.cos(z_rot), -math.sin(z_rot), 0],
                    [math.sin(z_rot),  math.cos(z_rot), 0],
                    [0, 0, 1]])
    

    R = R_x @ R_y @ R_x

    # transform coordinates
    node_coords = node_coords @ R.T

    # set new node coorinates
    for n, tag in enumerate(node_tags):
        gmsh.model.mesh.setNode(tag,node_coords[n],_)

    gmsh.model.occ.synchronize()


def isoScaleMesh(S):
    # ------------------------------------------------------------------------
    # Isotropic scaling of entire mesh
    # -------------------------------------------------------------------------

    gmsh.model.occ.synchronize()

    # get cartesian coordinates of all nodes    
    nodeTags, nodeCoord,_ = gmsh.model.mesh.getNodes()
    nodeCoords            = np.array(nodeCoord).reshape(-1,3)
    
    # transform coordinates
    nodeCoords = nodeCoords * S

    # set new node coorinates
    for n, nodeTag in enumerate(nodeTags):
        gmsh.model.mesh.setNode(nodeTag,nodeCoords[n],_)

    gmsh.model.occ.synchronize()
