"""
Geometric methods to compute volumes, areas, distances of the mesh entities
"""
# Geoemtric Module
# Create by Artur Castiel and Renata Tavares

#import pdb
import numpy as np
from numba import jit

# class for volume related problem
    #
    # def init_normal(self):
    #     self.core.create_tag_handle('NORMAL', 3)
    #     normal = np.zeros((len(self.core.all_faces), 3)).astype('float')
    #     index = 0
    #     for face in self.core.all_faces:
    #         verts = self.core.mb.get_connectivity(face)
    #         coords = np.array([self.core.mb.get_coords([vert]) for vert in verts])
    #         vec1 = coords[1] - coords[0]
    #         vec2 = coords[2] - coords[0]
    #         cross = np.cross(vec1,vec2)
    #         normal[index] = cross/np.linalg.norm(cross)
    #         index += 1
    #     self.core.set_data("NORMAL", normal, range_el=self.core.all_faces)
#@jit(parallel = True)

#@jit(parallel = True)
def normal_vec_2d(coords0,coords1):
    vec = coords1 - coords0
    norm = np.linalg.norm(vec, axis = 1)
    norm = 1/norm
    return np.array([vec[:,1], -vec[:,0], vec[:,2] ]).T * norm[:,np.newaxis]
    # distance = (np.inner(vec, vec, axis = 0))


# @jit(nopython=True)
# def cross_numba(vec1,vec2):
#     vec1 = double(vec1)
#     vec2 = double(vec2)
#     result = np.zeros((vec1.shape[0],3))
#     result[:,0] = vec1[:,1]*vec2[:,2] - vec1[:,2]*vec2[:,1]
#     result[:,1] = vec1[:,2]*vec2[:,0] - vec1[:,0]*vec2[:,2]
#     result[:,2] = vec1[:,0]*vec2[:,1] - vec1[:,1]*vec2[:,0]
#     return result
# @jit
# def cross_numba(vec1, vec2):
#     """ Calculate the cross product of two 3d vectors. """
#     result = np.zeros((vec1.shape[0],3)
#     # [a1, a2, a3] = double(vec1[0]), double(vec1[1]), double(vec1[2])
#     # [b1, b2, b3] = double(vec2[0]), double(vec2[1]), double(vec2[2])
#     # result[0] = a2 * b3 - a3 * b2
#     # result[1] = a3 * b1 - a1 * b3
#     # result[2] = a1 * b2 - a2 * b1
# return result
#@jit

def normal_vec(coords1, coords2, coords3):
    vec1 = coords1 - coords3
    vec2 = coords2 - coords3
    cross_product = np.cross(vec1, vec2)
    norm_cross = np.power(np.linalg.norm(cross_product,axis=1),-1)
    cross_product[:, 0] = norm_cross * cross_product[:, 0]
    cross_product[:, 1] = norm_cross * cross_product[:, 1]
    cross_product[:, 2] = norm_cross * cross_product[:, 2]
    return cross_product

#
# def normal_vec(coords0, coords1, coords2):
#     vec1 = coords1 - coords0
#     vec2 = coords2 - coords0
#     #cross = cross_numba(vec1,vec2)
#     cross = np.cross(vec1,vec2)
#     norm = np.linalg.norm(cross, axis = 1)
#     norm = 1/norm
#     return  cross * norm[:,np.newaxis]
#     # a = cross * norm

def point_distance(coords_1, coords_2):
    dist_vector = coords_1 - coords_2
    distance = sqrt(np.dot(dist_vector, dist_vector))
    return distance

#@jit(parallel=True)
#@jit
def get_average(coords_list):
    N = len(coords_list)
    return sum(coords_list)*(1/N)
#
# def tetraVolume(tet_nodes):
#     #input:
#     # A Matrix with 4x3 elements in which
#     # each line is one of the 4 nodes that
#     # a given tetrahedron is comprised
#     #ouput:
#     # the volume of the given tetrahedron
#     vect_1 = tet_nodes[1] - tet_nodes[0]
#     vect_2 = tet_nodes[2] - tet_nodes[0]
#     vect_3 = tet_nodes[3] - tet_nodes[0]
#     vol_eval = abs(np.dot(np.cross(vect_1, vect_2), vect_3))/6
#     return vol_eval
#
# def piramidVolume(pi_nodes):
#     #     P5           P4 _____ P3
#     #     /\             |     |
#     #    /  \            |     |
#     #   /____\           |_____|
#     # P1/4  P2/3        P1     P2
#     #input:
#     # A Matrix with 4x3 elements in which
#     # each line is one of the 5 nodes that
#     # a given piramid is comprised
#     # The 3 first nodes lie on the plane of the base (coplanar points) and must be connected.
#     # The fourth node of the matrix must be the top point (P5).
#     # ouput: the volume of the given piramid
#     vect_1 = pi_nodes[1] - pi_nodes[0]
#     vect_2 = pi_nodes[2] - pi_nodes[0]
#     base_area = abs(np.dot(vect_1, vect_2))
#
#     vect_3 = pi_nodes[3] - pi_nodes[0]
#     normal_vect = np.cross(vect_1, vect_2)
#     piram_height = np.dot(vect_3, normal_vect)/(np.linalg.norm(normal_vect)
#
#     piram_vol = (1/3)*base_area*piram_height
#     return(piram_vol)
#
# def hexahedronVolume(pi_nodes):
#     #
#     #    ______   <- F2
#     #   /     /|
#     #  /_____/ |
#     #  |     | |
#     #  | F1  | /
#     #  |_____|/
#     #
#     #
#     # P1 _____P2     P5  ____ P6
#     #   |     |         |    |
#     #   | F1  |         | F2 |
#     #   |_____|         |____|
#     # P4      P3     P8       P7
#     #
#     # F1 - Front Face
#     # F2 - Back Face
#     #  NOTE:
#     #  The given hexahedron may be irregular
#     #  The sketch above describes the connectivities
#     #
#     #input:
#     # A Matrix with 8x3 elements in which
#     # each line is one of the 8 nodes that
#     # a given hexahedron
#     #ouput:
#     # the volume of the given piramid
#     print(pi_nodes)
# def teste():
#     print("Entrou")
#     pass
