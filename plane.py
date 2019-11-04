import numpy as np
from preprocessor.meshHandle.configTools.configClass import coarseningInit as cm
from preprocessor.msCoarseningLib.partitionTools import partitionManager as pm
# from .preprocessor.meshHandle.finescaleMesh import FineScaleMesh

def plane_check(coords_list, elements, normal_plane, point_on_plane):


    pass


config_object = cm(empty=True)
config_object.smart(file='semi.msh')
former = pm(M, config_object)
former.run()



coord_list = M.nodes.coords[:]
elements = M.volumes.connectivities[:]
normal_plane = former.partitioner.dual.faces.normal[0]
point_on_plane = former.partitioner.dual.faces.center[0]
