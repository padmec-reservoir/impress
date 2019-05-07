"""
Module for management of fine scale mesh
"""
import time
import pdb
import numpy as np
from . corePymoab import CoreMoab
from . meshComponents import MeshEntities
import hdmedians as hd
import numpy as np
import scipy.spatial
print('Dual Coarse Mesh Module initialized')


class DualCoarseMesh:
    def __init__(self, M):
        self.M = M
        self.find_primal_coarse_centers()
        self.find_interface_centers()
        self.find_coarse_edges()
        self.find_coarse_faces()
        self.find_coarse_volumes()

    def find_primal_coarse_centers(self):
        self.coarse_center = np.zeros(len(self.M.coarse.elements)).astype("uint64")
        index = 0
        for coarse_volume in self.M.coarse.elements:
            center_coord = np.array(hd.geomedian(coarse_volume.faces.center[coarse_volume.faces.boundary].T))
            center_coord = center_coord.reshape((1,3))
            distance = scipy.spatial.distance.cdist(coarse_volume.volumes.center[:],center_coord)
            tag = np.nonzero(distance == distance.min())[0]
            if len(tag) != 1:
                tag = tag[0]
            center_volume = coarse_volume.volumes.all[tag]
            self.coarse_center[index] = coarse_volume.volumes.global_id[center_volume]
            index += 1

    def find_interface_centers(self):
        self.center_face = []
        coords = self.M.faces.center[:]
        for cface in self.M.coarse.interfaces_faces:
            tcoords = coords[cface].reshape(len(cface),3)
            u = int(cface[hd.medoid(tcoords,axis=0,indexonly=True)].ravel())
            #import pdb; pdb.set_trace()
            self.center_face.append(u)
            #self.center_face.append(cface[hd.medoid(tcoords,axis=0,indexonly=True) ] )

            #np.array(hd.geomedian(coarse_volume.faces.center[coarse_volume.faces.boundary].T))


        pass
    def find_coarse_edges(self):
        pass
    def find_coarse_faces(self):
        pass
    def find_coarse_volumes(self):
        pass

class GraphMesh:
    def __init__(self, M):
        print("Teste")
        print(M)
