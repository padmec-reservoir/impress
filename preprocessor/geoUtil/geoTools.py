"""
Geometric methods to compute volumes, areas, distances of the mesh entities
"""
# Geoemtric Module
# Created by Artur Castiel, Renata Tavares and Filipe Cumaru.

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

def cross_product_3d(u, v):
    w = np.zeros(3)
    w[0] = u[1]*v[2] - u[2]*v[1]
    w[1] = u[2]*v[0] - u[0]*v[2]
    w[2] = u[0]*v[1] - u[1]*v[0]
    return w

def triangle_area(v1, v2, v3):
    w = cross_product_3d(v1 - v2, v1 - v3)
    area = 0.5*np.linalg.norm(w)
    return area

def polygon_area(moab_core, polygon):
    """
    Calculate the area of a polygon by triangulation.
    """
    # Retrieve vertices handles and coordinates from face handle.
    vertices = moab_core.get_adjacencies(polygon, 0)
    vert_coords = moab_core.get_coords(vertices).reshape(len(vertices), 3)
    vertices_dict = dict(zip(vertices, vert_coords))

    # If the polygon is a triangle, then just compute the area by
    # definition.
    if moab_core.type_from_handle(polygon) == types.MBTRI or vertices.size() == 3:
        return triangle_area(vert_coords[0], vert_coords[1], vert_coords[2])
    
    # Else, compute a triangulation for this shape.
    mtu = topo_util.MeshTopoUtil(moab_core)
    # Choose a vertex to start and compute its neighbors, a.k.a, the vertices
    # sharing an edge with it.
    v0 = vertices[0]
    v0_neighbors = rng.intersect(vertices, mtu.get_bridge_adjacencies(v0, 1, 0))
    v0_coords = vertices_dict[v0]

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
        vj = rng.subtract(vi_neighbors, visited_verts)[0]
        vi_coords, vj_coords = vertices_dict[vi], vertices_dict[vj]
        area += triangle_area(v0_coords, vi_coords, vj_coords)
        visited_verts.insert(vj)
        vi = vj
    
    return area
    
def pyramid_volume(moab_core, face, v):
    """
    Compute the volume of a pyramid.

    moab_core: a PyMOAB core instance.
    face: a PyMOAB EntityHandle representing the base of the pyramid.
    v: a NumPy array representing the coordinates of the top vertex.
    """
    # Get three vertices from the base.
    vertices_handles = moab_core.get_connectivity(face)[0:3]
    p1, p2, p3 = moab_core.get_coords(vertices_handles).reshape((3,3))

    # Compute the area of the base.
    A = polygon_area(moab_core, face)

    # Compute the distance from v to the base plane.
    u = cross_product_3d(p2 - p1, p3 - p1)
    n = u / np.linalg.norm(u)
    h = np.abs(n.dot(v - p1))

    return (A*h) / 3

def polyhedron_volume(moab_core, polyhedron, center):
    """
    Computes the volume of a convex polyhedron.

    moab_core: a PyMOAB core instance.
    face: a PyMOAB EntityHandle representing the polyhedron.
    center: a NumPy array containing the coordinates of the
            centroid.
    """
    faces = moab_core.get_adjacencies(polyhedron, 2)
    vertices = moab_core.get_adjacencies(polyhedron, 0)

    # If the polyhedron is a pyramid or a tetrahedron, then compute
    # its volume straight ahead.
    if moab_core.type_from_handle(polyhedron) == types.MBTET:
        base = faces[0]
        base_vertices = moab_core.get_adjacencies(base, 0)
        top_vertex = rng.subtract(vertices, base_vertices)[0]
        top_vertex_coords = moab_core.get_coords(top_vertex)
        volume = pyramid_volume(moab_core, base, top_vertex_coords)
    elif moab_core.type_from_handle(polyhedron) == types.MBPYRAMID:
        base = [face for face in faces if moab_core.type_from_handle(face) != types.MBTRI][0]
        base_vertices = moab_core.get_adjacencies(base, 0)
        top_vertex = rng.subtract(vertices, base_vertices)[0]
        top_vertex_coords = moab_core.get_coords(top_vertex)
        volume = pyramid_volume(moab_core, base, top_vertex_coords)
    # Otherwise, compute the volume by splitting the polyhedron
    # into pyramids, each face acting like the base and the centroid
    # acting like the top vertex.
    else:
        volume = sum([pyramid_volume(moab_core, face, center) for face in faces])

    return volume
