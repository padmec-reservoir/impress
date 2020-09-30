"""
Geometric methods to compute volumes, areas, distances of the mesh entities
"""
# Geoemtric Module
# Created by Artur Castiel and Renata Tavares

import numpy as np
from pymoab import topo_util, types, rng

def normal_vec_2d(coords0, coords1):
    vec = coords1 - coords0
    vec_norm = np.linalg.norm(vec, axis = 1)
    return np.array([vec[:,1], -vec[:,0], vec[:,2] ]).T / vec_norm[:,np.newaxis]

def normal_vec(coords1, coords2, coords3):
    vec1 = coords1 - coords3
    vec2 = coords2 - coords3
    cross_product = np.cross(vec1, vec2)
    norm_cross = np.power(np.linalg.norm(cross_product,axis=1),-1)
    cross_product[:, 0] = norm_cross * cross_product[:, 0]
    cross_product[:, 1] = norm_cross * cross_product[:, 1]
    cross_product[:, 2] = norm_cross * cross_product[:, 2]
    return cross_product

def point_distance(coords_1, coords_2):
    dist_vector = coords_1 - coords_2
    distance = np.sqrt(np.dot(dist_vector, dist_vector))
    return distance

def get_average(coords_list):
    N = len(coords_list)
    return sum(coords_list)*(1/N)

def triangle_area(v1, v2, v3):
    u = np.cross(v1 - v2, v1 - v3)
    area = 0.5*np.linalg.norm(u)
    return area

def polygon_area(moab_core, polygon):
    """
    Calculate the area of a polygon by triangulation.
    """
    # Retrieve vertices handles and coordinates from face handle.
    vertices = rng.Range(moab_core.get_connectivity(polygon))

    # If the polygon is a triangle, then just compute the area by
    # definition.
    if moab_core.type_from_handle(polygon) == types.MBTRI or vertices.size() == 3:
        vert_coords = moab_core.get_coords(vertices)
        vert_coords = vert_coords.reshape(3, 3)
        return triangle_area(vert_coords[0], vert_coords[1], vert_coords[2])
    
    # Else, compute a triangulation for this shape.
    mtu = topo_util.MeshTopoUtil(moab_core)
    # Choose a vertex to start and compute its neighbors, a.k.a, the vertices
    # sharing an edge with it.
    v0 = vertices[0]
    v0_neighbors = rng.intersect(vertices, mtu.get_bridge_adjacencies(v0, 1, 0))
    v0_coords = moab_core.get_coords(v0)

    # vi is the vertice currently being visited, and vj is its neighbor not
    # visited yet. At each iteration, we compute the area of the triangle (v0, vi, vj).
    vi, vj = v0_neighbors[0], None

    # At each iteration, we store the vertices that already took part in 
    # a triangle to avoid recalculating the same triangle.
    visited_verts = rng.Range([v0, vi])

    # While there still vertices to be visited, compute a new triangle.
    area = 0
    while visited_verts.size() < vertices.size():
        vi_neighbors = rng.intersect(vertices, mtu.get_bridge_adjacencies(vi, 1, 0))
        vj = [v for v in vi_neighbors if v not in visited_verts][0]
        vi_coords, vj_coords = moab_core.get_coords([vi, vj]).reshape((2,3))
        area += triangle_area(v0_coords, vi_coords, vj_coords)
        visited_verts.insert(vj)
        vi = vj
    
    return area
    

def polyhedron_volume(moab_core, polyhedron):
    pass
