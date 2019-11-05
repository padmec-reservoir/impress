import numpy as np
from preprocessor.meshHandle.configTools.configClass import coarseningInit as cm
from preprocessor.msCoarseningLib.partitionTools import partitionManager as pm
# from .preprocessor.meshHandle.finescaleMesh import FineScaleMesh
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
import preprocessor.meshHandle.imprutil as ipw
from scipy.spatial import ConvexHull
slen = np.vectorize(len)


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
config_object.smart(file='semi.msh')
former = pm(M, config_object)
former.run()

x = 70
coord_list = M.nodes.coords[:]
elements = M.volumes.connectivities[:]
normal_plane = former.partitioner.dual.faces.normal[x]


point_on_plane = former.partitioner.dual.faces.center[x]
nodes_faces = former.partitioner.dual.faces.connectivities[x]
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
