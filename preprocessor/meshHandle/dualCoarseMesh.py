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
import scipy as sp
import scipy.spatial
import collections

print('Dual Coarse Mesh Module initialized')


class DualCoarseMesh:
    def __init__(self, M):
        self.M = M
        self.find_interface_centers()
        self.find_primal_coarse_centers()
        self.find_vol_neighbors_to_interface_center()
        # pdb.set_trace()

        self.find_coarse_edges()
        self.find_coarse_faces()
        # self.find_coarse_volumes()

    def find_primal_coarse_centers(self):
        self.coarse_center = np.zeros(len(self.M.coarse.elements)).astype("uint64")
        index = 0
        for index, coarse_volume in enumerate(self.M.coarse.elements):
            center_coord = np.array(hd.geomedian(coarse_volume.faces.center[coarse_volume.faces.boundary].T))
            center_coord = center_coord.reshape((1,3))
            distance = scipy.spatial.distance.cdist(coarse_volume.volumes.center[:],center_coord)
            tag = np.nonzero(distance == distance.min())[0]
            if len(tag) != 1:
                tag = tag[0]
            center_volume = coarse_volume.volumes.all[tag]
            self.coarse_center[index] = coarse_volume.volumes.father_id[center_volume]
            #index += 1

    def find_interface_centers(self):
        #self.center_face = []
        coords = self.M.faces.center[:]
        self.center_face = np.zeros((len(self.M.coarse.interfaces_faces), 1)).astype("uint64")
        index = 0
        for cface in self.M.coarse.interfaces_faces:
            tcoords = coords[cface].reshape(len(cface),3)
            center_coord = np.array(hd.geomedian(tcoords.T)).reshape((1,3))
            distance = scipy.spatial.distance.cdist(tcoords, center_coord)
            tag = np.nonzero(distance == distance.min())[0]
            if len(tag) != 1:
                tag = tag[0]
            self.center_face[index] = int(cface[tag].ravel())
            index += 1


    def find_vol_neighbors_to_interface_center(self):
        # internal faces
        internal_volumes = self.M.faces.bridge_adjacencies(self.center_face[0:self.M.coarse.num_internal_faces], interface ="faces", target = "volumes")
        external_volumes = self.M.faces.bridge_adjacencies(self.center_face[self.M.coarse.num_internal_faces:], interface ="faces", target = "volumes")
        self.interface_vol = np.vstack((internal_volumes,np.hstack((external_volumes,external_volumes))))
        tag = self.M.coarse.partition[:].ravel()[self.interface_vol][:,0] > self.M.coarse.partition[:].ravel()[self.interface_vol][:,1]
        #pdb.set_trace()
        self.interface_vol[tag,0], self.interface_vol[tag,1] = self.interface_vol[tag,1], self.interface_vol[tag,0]
        #pdb.set_trace()

    def find_coarse_edges(self):
        edges = np.array([])
        coarse_edges = [np.array([])]* len(self.M.coarse._faces)
        for x in range(len(self.M.coarse.elements)):
            cneigh, cfaces  = self.M.coarse.iface_neighbors(x)
            pointer = self.M.coarse.partition[:].ravel()[self.interface_vol[cfaces,:]] == x
            pointer_to_face = cneigh >= len(self.M.coarse.elements)
            pointer[pointer_to_face, :] = False
            target_volume = self.interface_vol[cfaces,:][pointer]
            target_volume_inside = self.interface_vol[cfaces,:][pointer_to_face,0]
            element_target = self.M.coarse.father_to_local_id(target_volume, "volumes", x).ravel()
            element_target_inside = self.M.coarse.father_to_local_id(target_volume_inside, "volumes", x).ravel()
            element_center = self.M.coarse.father_to_local_id(self.coarse_center[x], "volumes", x)
            shortest = GraphMesh(self.M.coarse.elements[x],center = element_center)
            # if x == 22:
            #     pdb.set_trace()
            for index, el in enumerate(element_target):

                # if index == 3 and x == 22:
                #     pdb.set_trace()
                #     print(index)
                center = self.M.coarse.elements[x].volumes.father_id[el]
                face = np.logical_or((self.interface_vol[:, 0] == center), (self.interface_vol[:, 1] == center))
                face = (np.where(face)[0])
                line = self.interface_vol[face,:]
                tag = np.zeros(line.shape[0],dtype=bool)
                #pdb.set_trace()
                tag = line == center
                for index, f in enumerate(face):
                    t = np.where(tag[index,:])[0]
                    m = self.M.coarse.elements[x].volumes.father_id[shortest.path(el)]
                    if t == 0:
                        coarse_edges[f] = np.append(m, coarse_edges[f])
                    else:
                        coarse_edges[f] = np.append(coarse_edges[f], m)

            for index, el in enumerate(element_target_inside):
                center = self.M.coarse.elements[x].volumes.father_id[el]
                face = np.logical_or((self.interface_vol[:, 0] == center), (self.interface_vol[:, 1] == center))
                face = int(np.where(face)[0])
                m = self.M.coarse.elements[x].volumes.father_id[shortest.path(el)]
                coarse_edges[face] = np.append(m, coarse_edges[face])
            #pdb.set_trace()
            self.coarse_edges = np.array(coarse_edges)

    def find_coarse_faces(self):
        pass
    def find_coarse_volumes(self):
        pass


class GraphMesh:
    def __init__(self, N, center):
        graph = self.create_sparse_matrix(N.faces.bridge_adjacencies(N.faces.internal, interface="faces", target="volumes"), len(N.volumes))
        # self.dist_matrix, self.predecessors = sp.sparse.csgraph.shortest_path(graph, directed=False, return_predecessors = True)
        self.dist_matrix, self.predecessors = sp.sparse.csgraph.dijkstra(graph, directed=False, return_predecessors = True, indices = center)
        self.predecessors = self.predecessors.ravel()

#         self.indicies = target
#         self.center = center
#         self.cedges = []
#         for el in target:
#             self.cedges.append(self.path(el).ravel())
# #             np.concatenate((self.cedges, self.path(el).ravel()))

    def create_sparse_matrix(self, graph_edges, size_vol):
        sparse_matrix = sp.sparse.coo_matrix((np.ones((len(graph_edges),),dtype=bool).T, (graph_edges[:, 0], graph_edges[:, 1])), shape=(size_vol, size_vol))
        return sparse_matrix

    def path(self, target):
        index = target
        path = []
        while index != -9999:
            path.append(index)
            index = self.predecessors[index]
        return np.array(path)
