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
        self.find_primal_coarse_centers()
        self.find_interface_centers()
        self.find_vol_neighbors_to_interface_center()
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



        # self.interface_vol = (np.hstack((internal_volumes, external_volumes)).astype("uint64"))
        # face_classification = np.zeros((len(internal_volumes) ,2))
        # face_classification[:,0] , face_classification[:,1] = np.arange(0,len(internal_volumes)), np.arange(0,len(internal_volumes))
        # face_classification = np.hstack((face_classification.ravel(), np.arange(len(internal_volumes), len(internal_volumes) + len(external_volumes))))
        # self.interface_vol, indices = np.unique(self.interface_vol, return_index=True)
        # self.face_classification = face_classification[indices].astype("uint64")


    def find_coarse_edges(self):
        cneigh = self.M.coarse.iface_neighbors
        coarse_reference = np.zeros(cneigh.shape)
        part = self.M.coarse.partition[:]

        pdb.set_trace()
        coarse_reference[:,0], coarse_reference[:,1] = self.M.coarse.partition[self.interface_vol[:,0]].ravel(), self.M.coarse.partition[self.interface_vol[:,1]].ravel()
        edges = np.array([])

        coarse_edges = [collections.deque() for el in range(len(self.M.coarse._faces))]



        for x in range(len(self.M.coarse.elements)):
            element_target= self.M.coarse.father_to_local_id(self.interface_vol[np.where(coarse_reference == x)[0]], "volumes", x).ravel()
            face_target = self.face_classification[np.where(coarse_reference == x)[0]].astype("uint64")
            element_center = self.M.coarse.father_to_local_id(self.coarse_center[x], "volumes", x)
            pdb.set_trace()
            faces = self.M.coarse.faces_neighbors[x,np.where(self.M.coarse.faces_neighbors[x,:] != None)[0]].astype("uint64")



            shortest = GraphMesh(self.M.coarse.elements[x], target = element_target, center = element_center)
            edges = np.concatenate((edges,self.M.coarse.elements[x].volumes.father_id[shortest.cedges]))


            #[edges.append(shortest.path(element_center, el)) for el in element_target]

        self.coarse_edges = np.unique(edges).astype("uint64")

    def find_coarse_faces(self):
        pass
    def find_coarse_volumes(self):
        pass


class GraphMesh:
    def __init__(self, N, target, center):
        graph = self.create_sparse_matrix(N.faces.bridge_adjacencies(N.faces.internal, interface="faces", target="volumes"), len(N.volumes))
        # self.dist_matrix, self.predecessors = sp.sparse.csgraph.shortest_path(graph, directed=False, return_predecessors = True)
        self.dist_matrix, self.predecessors = sp.sparse.csgraph.dijkstra(graph, directed=False, return_predecessors = True, indices = center)
        self.predecessors = self.predecessors.ravel()
        self.indicies = target
        self.center = center
        self.cedges = np.array([])
        for el in target:
            self.cedges = np.concatenate((self.cedges, self.path(el).ravel()))
        self.cedges = np.unique(self.cedges)

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
