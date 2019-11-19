import numpy as np
from preprocessor.meshHandle.configTools.configClass import coarseningInit as cm
from preprocessor.msCoarseningLib.partitionTools import partitionManager as pm
# from .preprocessor.meshHandle.finescaleMesh import FineScaleMesh
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
import preprocessor.meshHandle.imprutil as ipw
from scipy.spatial import ConvexHull
slen = np.vectorize(len)

def global_to_local(global_matrix):
    global_to_local_dic = {}
    global_connectivities_id = np.unique(global_matrix)
    local_connectivities_id = np.arange(len(global_connectivities_id))
    for key, val in zip(global_connectivities_id,
                        local_connectivities_id):
        global_to_local_dic[key] = val
    return np.vectorize(global_to_local_dic.get)(global_matrix)


def check_nodes_in_volume(el_coords, faces_connectivities, faces_normal, centers):
    el_center = el_coords.mean(axis=0).reshape((-1, 3))
    import pdb; pdb.set_trace()
    faces_center = el_coords[faces_connectivities].mean(axis=1)
    pseudo_normal = faces_center - el_center.repeat(len(faces_center), axis=0)
    change_sign = (faces_normal*pseudo_normal).sum(axis=1) <= 0

    faces_normal[change_sign] = -1 * faces_normal[change_sign]
    # o erro esta aqui

    nodes_indicator = np.zeros((len(centers),
                                faces_connectivities.shape[0]), dtype=bool)
    for el in range(faces_connectivities.shape[0]):
        plane_check = semi_plan_check(centers, faces_normal[el],
                             el_coords[faces_connectivities[el,0]], tol=1e-8)
        nodes_indicator[:, el] = plane_check
    import pdb; pdb.set_trace()
    return nodes_indicator.any(axis=1)


def semi_plan_check(coords_list, normal_plane, point_on_plane, tol=1e-20):
    """coord_list -> list of points to be checked
    normal_plane -> normal vector to plane used for the comparision
    point_on_plane -> a point on the plane
    returns a bool of the size of coords_list with True or false """
    center_to_coords = coords_list - \
        np.repeat(point_on_plane.reshape((-1,3)), len(coords_list), axis=0)
    normal_plane = \
        np.repeat(normal_plane.reshape((-1,3)), len(coords_list), axis=0)
    inner_product = np.sum(center_to_coords*normal_plane,axis=1)
    flag = np.zeros(inner_product.shape, dtype=bool)
    flag[inner_product >= 0] = True

    #flag[np.abs(inner_product) < tol] =True
    return flag


def plane_check(coords_list, elements, normal_plane, point_on_plane, tol=1e-20):
    # import pdb; pdb.set_trace()
    center_to_coords = coord_list - \
        np.repeat(point_on_plane.reshape((-1,3)), len(coords_list), axis=0)
    normal_plane = \
        np.repeat(normal_plane.reshape((-1,3)), len(coords_list), axis=0)
    inner_product = np.sum(center_to_coords*normal_plane,axis=1)
    flag = np.zeros(inner_product.shape, dtype=int)
    flag[inner_product > 0] = 1
    flag[inner_product < 0] = -1
    flag[np.abs(inner_product) < tol] = 0
    # qq = np.array([np.array([1,3,4]),np.array([1,2])])
    index = elements.shape[1]
    flag = np.abs(np.sum(flag[elements], axis=1))
    cross_elements = np.zeros((len(elements), 1), dtype=bool)
    cross_elements[flag < (index-1)] = True
    return cross_elements.T.ravel()

#M = msh('mesh/semi3.msh', dim = 3)

config_object = cm(empty=True)
config_object.smart(file='semi2.msh')
former = pm(M, config_object)
former.run()

x = 70
coord_list = M.nodes.coords[:]
elements = M.volumes.connectivities[:]
normal_plane = former.partitioner.dual.faces.normal[x]


point_on_plane = former.partitioner.dual.faces.center[x]
nodes_faces = former.partitioner.dual.faces.connectivities[x]


el_coord = former.partitioner.primal.nodes.coords[former.partitioner.primal.volumes.connectivities[x]]
elements_coord = global_to_local(former.partitioner.primal.volumes.connectivities[x])
lag = former.partitioner.primal.volumes.adjacencies[x]
faces_connectivities = global_to_local(former.partitioner.primal.faces.connectivities[former.partitioner.primal.volumes.adjacencies[x]])

faces_normal = np.zeros((len(lag), 3))
# norm_tmp = np.vectorize(former.partitioner.primal.faces._normal)
# faces_normal = norm_tmp(lag)
for el in range(len(faces_normal)):
    faces_normal[el, :] = former.partitioner.primal.faces._normal(lag[el])


lemo = check_nodes_in_volume(el_coord, faces_connectivities, faces_normal, M.volumes.center[:])


#import pdb; pdb.set_trace()
qq = semi_plan_check(M.volumes.center[:], normal_plane, point_on_plane, tol=1e-20)
qq = plane_check(coord_list, elements, normal_plane, point_on_plane)

normal_plane_rep = np.repeat(normal_plane, len(nodes_faces), axis=0)
gam = 1.5
top_coord = coord_list[nodes_faces] + (normal_plane_rep*gam)
bot_coord = coord_list[nodes_faces] - (normal_plane_rep*gam)

all_coords = np.vstack((top_coord, bot_coord))
all_coords = M.nodes.coords[M.nodes.boundary]

hull = ConvexHull(all_coords)
vec_test = ipw.point_in_volumes(all_coords, hull.simplices.astype('int64') , M.volumes.center[qq], 0)


dd = np.where(qq)[0][vec_test!=0]

M.pressure[:] = 0
M.pressure[qq] = 1
# M.pressure[dd] = 2
M.core.print(file='plane',case='test')
