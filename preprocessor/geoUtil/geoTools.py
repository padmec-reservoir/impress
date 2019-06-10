"""
Geometric methods to compute volumes, areas, distances of the mesh entities
"""
# Geoemtric Module
# Create by Artur Castiel and Renata Tavares

import numpy as np
from numpy import linalg as la
from numba import jit
#import pdb

# Calculate areas

def triangle_area(p0, p1, p2):
    """
         P0
        / \
       /   \
      /_ _ _\
     P1      P2

Calculate the area of a single triangle or a group of triangles giving the coordinates of its points.

-> Input:
Pn for n = 0,1,2 are arrays containing the coordinates of a given point of a convex triangle. Pn has length mx3, where m is the number of triangles given. The order of the points is not important in this function

-> Output
A row vector containing the area of the given triangles.                """

    # Setting the vectors
    v1 = p1 - p0
    v2 = p2 - p0
    tri_area = (1/2)*la.norm(np.cross(v1, v2), axis = 1)
    return tri_area

def quadrilateral_area(p0, p1, p2, p3):
    """
     P0 ____ P2
       |    |
       |    |
       |____|
     P1      P3

Calculate the area of a single quadrilateral or a group of quadrilaterals giving the coordinates of its points.

-> Input:
Pn for n = 1,2...4 are arrays containing the coordinates of a given point of a quadrilateral. Pn has length mx3, where m is the number of quadrilaterals given. The sequence of the nodes in the matrix must be the same as the figure.

-> Output:
A row vector containing the volume of the volumes of the given hexahedrons """

    # Setting the diagonals
    d1 = p0 - p3
    d2 = p1 - p2
    quad_area = (1/2)*la.norm(np.cross(d1, d2), axis = 1)
    return quad_area

# Calculate volumes

def pentagon_area(p0, p1, p2, p3, p4):
    # Implementation of Sotke's theorem
    # The order of the point is important

    """
      P0
      /\
     /  \
 P1 /    \ P4
    \    /
     \__/
    P2  P3

    """
    coords = np.array([p0, p1, p2, p3, p4])
    sum = [0,0,0]
    for i in range(len(coords)):
        c1 = poly[i]
        c2 = poly[(i+1) % N]
        prod = np.cross(coords[i], [(i+1) % N])
        sum[0] += prod[0]
        sum[1] += prod[1]
        sum[2] += prod[2]
    pentagon_area = (1/2)*np.dot(total, normal_vec(p0, p1, p2))
    return pentagon_area

def hexagon_area(p0, p1, p2, p3, p4, p5):
    # Implementation of Sotke's theorem
    """
          _____
         /     \
        /       \
        \       /
         \_____/

    """
    coords = np.array([p0, p1, p2, p3, p4, p5])
    sum = [0,0,0]
    for i in range(len(coords)):
        c1 = poly[i]
        c2 = poly[(i+1) % N]
        prod = np.cross(coords[i], [(i+1) % N])
        sum[0] += prod[0]
        sum[1] += prod[1]
        sum[2] += prod[2]
    hexagon_area = (1/2)*np.dot(total, normal_vec(p0, p1, p2))
    return hexagon_area

def square_piramid_volumes(p0, p1, p2, p3, p4):
    """
   P1 ______ P3
     |\    /|
     | \  / |
     |  \/  |
     |  /\P0|
     | /  \ |
     |/____\|
   P2        P4

Calculate the volume of a single pyramid or a group of pyramids (its base must be a quadrilateral) giving the coordinates of its points.

-> Input:
Pn for n = 0, 1, 2... 4 are arrays containing the coordinates of a given point of a pymarid. Pn has length mx3, where m is the number of pyramids given. The sequence of the nodes in the matrix must be the same as the figure.

-> Output:
A row vector containing the volume of the volumes of the given pyramids"""

    base_area = quadrilateral_area(p1, p2, p3, p4)
    v1 = p1 - p0
    v2 = p1 - p2
    v3 = p1 - p3
    h = np.linalg.norm(np.linalg.dot(v1,np.linalg.cross(v2,v2))/np.linalg.norm(np.linalg.cross(v2,v2))
    piramid_vol = (1/3)*base_area*h
    return piramid_vol

def pentagon_piramid_volumes(p0, p1, p2, p3, p4, p5):
    """
DRAW

Calculate the volume of a single pyramid or a group of pyramids (its base must be a quadrilateral) giving the coordinates of its points.

-> Input:
Pn for n = 0, 1, 2... 4 are arrays containing the coordinates of a given point of a pymarid. Pn has length mx3, where m is the number of pyramids given. The sequence of the nodes in the matrix must be the same as the figure.

-> Output:
A row vector containing the volume of the volumes of the given pyramids"""

    base_area = pentagon_area(p1, p2, p3, p4, p5)
    v1 = p1 - p0
    v2 = p1 - p2
    v3 = p1 - p3
    h = np.linalg.norm(np.linalg.dot(v1,np.linalg.cross(v2,v2))/np.linalg.norm(np.linalg.cross(v2,v2))
    piramid_vol = (1/3)*base_area*h
    return piramid_vol

def hexagon_piramid_volumes(p0, p1, p2, p3, p4, p5, p6):
    """

DRAW

Calculate the volume of a single pyramid or a group of pyramids (its base must be a quadrilateral) giving the coordinates of its points.

-> Input:
Pn for n = 0, 1, 2... 4 are arrays containing the coordinates of a given point of a pymarid. Pn has length mx3, where m is the number of pyramids given. The sequence of the nodes in the matrix must be the same as the figure.

-> Output:
A row vector containing the volume of the volumes of the given pyramids"""

    base_area = hexagon_area(p1, p2, p3, p4, p5, p6)
    v1 = p1 - p0
    v2 = p1 - p2
    v3 = p1 - p3
    h = np.linalg.norm(np.linalg.dot(v1,np.linalg.cross(v2,v2))/np.linalg.norm(np.linalg.cross(v2,v2))
    piramid_vol = (1/3)*base_area*h
    return piramid_vol


def tetrahedron_volume(p0, p1,p2,p3): # Documentar
    """
          P0
          /|\
         / | \
        /  |  \
   P1 / _ _|_ _ \  P2
      \    |    /
        \  |  /
          \|/
           P3

Calculate the volume of a single tetrahedron or a group of tetrahedrons giving the coordinates of its points.

-> Input:
Pn for n = 1,2...4 are arrays containing the coordinates of a given point of a quadrilateral. Pn has length mx3, where m is the number of tetrahedrons given. The order of the points is not important in this function.

-> Output:
A row vector containing the volume of the volumes of the given tetrahedrons """
    # Setting the vectors
    v1 = p0 - p1
    v2 = p0 - p2
    v3 = p0 - p3
    tetra_vol = (1/6)*np.sum((np.cross(v1,v2))*v3, axis = 1)
    return tetra_vol

def hexahedron_volume(p0, p1, p2, p3, p4, p5, p6, p7):
    """ ______
       /     /|
      /_____/ |
      |     | |    <- F2
      | F1  | /
      |_____|/


     P4 ____ P5     P6 ____ P7
       |    |         |    |
       | F1 |         | F2 |
       |____|         |____|
     P0      P1     P2       P3

     F1 - Front Face
     F2 - Back Face
      NOTE:
     The given hexahedron may be irregular
     The method used here includes this possibility
     The sketch above describes the connectivities

Calculate the volume of a single hexahedron or a group of hexahedrons giving the coordinates of its points.

-> Input:
Pn for n = 1,2...8 are arrays containing the coordinates of a given point of a hexahedron. Pn has length mx3, where m is the number of hexahedrons given. The sequence of the nodes in the matrix must be the same as the figure.

-> Ouput:
A row vector containing the volume of the volumes of the given hexahedrons                                                 """

    v1 = np.cross((p7-p0), (p1-p0))
    v2 = np.cross((p7-p0), (p4-p0))
    v3 = np.cross((p7-p0), (p2-p0))

    vol1 = np.sum(v1*(p3-p5), axis = 1)
    vol2 = np.sum(v2*(p5-p6), axis = 1)
    vol3 = np.sum(v3*(p6-p3), axis = 1)

    hexa_vol = (vol1+vol2+vol3)/6
    return hexa_vol

# Calculate things

@jit
def normal_vec(coords0, coords1, coords2):
    vec1 = coords1 - coords0
    vec2 = coords2 - coords0
    #cross = cross_numba(vec1,vec2)
    cross = np.cross(vec1,vec2)
    norm = np.linalg.norm(cross, axis = 1)
    norm = 1/norm
    return  cross * norm[:,np.newaxis]
    # a = cross * norm

def point_distance(coords_1, coords_2):
    dist_vector = coords_1 - coords_2
    distance = sqrt(np.dot(dist_vector, dist_vector))
    return distance

#@jit(parallel=True)
@jit
def get_average(coords_list):
    N = len(coords_list)
    return sum(coords_list)*(1/N)

@jit(parallel = True)
def normal_vec_2d(coords0,coords1):
    vec = coords1 - coords0
    norm = np.linalg.norm(vec, axis = 1)
    norm = 1/norm
    return np.array([vec[:,1], -vec[:,0], vec[:,2] ]).T * norm[:,np.newaxis]
    # distance = (np.inner(vec, vec, axis = 0))

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
