import numpy as np
from preprocessor.meshHandle.configTools.configClass import coarseningInit as cm
from preprocessor.msCoarseningLib.partitionTools import partitionManager as pm
# from .preprocessor.meshHandle.finescaleMesh import FineScaleMesh

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
    cross_elements[flag != index] = True
    return cross_elements.T.ravel()

config_object = cm(empty=True)
config_object.smart(file='semi.msh')
former = pm(M, config_object)
former.run()

coord_list = M.nodes.coords[:]
elements = M.volumes.connectivities[:]
normal_plane = former.partitioner.dual.faces.normal[27]
point_on_plane = former.partitioner.dual.faces.center[27]

qq = plane_check(coord_list, elements, normal_plane, point_on_plane)
M.pressure[:] = 0
M.pressure[qq] = 1
M.core.print(file='plane',case='test')
