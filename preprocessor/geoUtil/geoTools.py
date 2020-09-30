"""
Geometric methods to compute volumes, areas, distances of the mesh entities
"""
# Geoemtric Module
# Created by Artur Castiel and Renata Tavares

import numpy as np

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
