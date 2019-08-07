"""
Use of Pymoab methods to manage the multiscale mesh
"""
#import pdb
from . corePymoab import CoreMoab
from pymoab import core, types, rng, topo_util
from pymoab import skinner as sk
import numpy as np
import time

class MsCoreMoab(CoreMoab):
    def __init__(self, father_core, num, coarse_vec):
        self.father_core = father_core
        self.dimension = father_core.dimension
        self.mb = father_core.mb
        self.level = father_core.level + 1
        self.coarse_num = num
        self.father_root_set = father_core.root_set
        self.root_set = self.mb.create_meshset(types.MESHSET_TRACK_OWNER) #types.MESHSET_TRACK_OWNER)
        self.mtu = father_core.mtu
        self.handleDic = father_core.handleDic
        if self.dimension == 3:
            self.all_volumes = self.range_index(coarse_vec, range_handle=father_core.all_volumes)
            self.all_faces = self.mb.get_adjacencies(self.all_volumes, 2, False, op_type=types.UNION)
            self.all_edges = self.mb.get_adjacencies(self.all_volumes, 1, False, op_type=types.UNION)
            self.all_nodes = self.mb.get_adjacencies(self.all_volumes, 0, False, op_type=types.UNION)
        elif self.dimension == 2:
            self.all_volumes = rng.Range()
            self.all_faces = self.range_index(coarse_vec, range_handle=father_core.all_faces)
            self.all_edges = self.access_handle(self.all_faces)
            self.all_nodes = rng.Range(self.mb.get_connectivity(self.all_faces))
        self.mb.add_entities(self.root_set, self.all_volumes)
        self.mb.add_entities(self.root_set, self.all_faces)
        self.mb.add_entities(self.root_set, self.all_edges)
        self.mb.add_entities(self.root_set, self.all_nodes)
        all_entities = self.mb.get_entities_by_handle(self.root_set)

        [self.boundary_nodes, self.boundary_edges, self.boundary_faces, self.boundary_volumes] = self.skinner_operation()
        self.internal_nodes = rng.subtract(self.all_nodes, self.boundary_nodes)
        self.internal_edges = rng.subtract(self.all_edges, self.boundary_edges)
        self.internal_faces = rng.subtract(self.all_faces, self.boundary_faces)
        self.internal_volumes = rng.subtract(self.all_volumes, self.boundary_volumes)


        # self.print_set2 = self.mb.create_meshset()
        # self.mb.add_entities(self.print_set2, self.all_faces)
        # self.mb.write_file('auxfaces'+f'{self.coarse_num}'+'.h5m', [self.print_set2])


        if self.level == 1:
            self.id_name = "LOCAL_ID_L" + str(self.level) + "-" + str(self.coarse_num)
            self.father_id_name = self.father_core.id_name
        else:
            self.father_id_name = self.father_core.id_name
            self.id_name = self.father_id_name + str("L") + str(self.level) + "-" + str(self.coarse_num)


        self.init_id()
        # all_entities = self.mb.get_entities_by_handle(self.root_set)
        # print(all_entities)
        self.flag_dic = {key:[rng.intersect(all_entities,el) for el in value] for (key, value) in father_core.flag_dic.items()}

    def auxskinner_operation(self):
        #skin = sk.Skinner(self.mb)
        # print("Entering skinner test")

        if self.dimension == 3:
            faces_on_skin_handles = self.bridge_adjacencies(self.all_faces,self.dimension)
            edges_on_skin_handles = self.access_handle(faces_on_skin_handles)
            nodes_on_skin_handles = self.access_handle(edges_on_skin_handles)
            nodes_in_volumes = ([self.mb.get_adjacencies(el_handle,0) for el_handle in self.all_volumes])
            check_volumes = ([(rng.intersect(el_handle,nodes_on_skin_handles))  for el_handle in nodes_in_volumes])
            external_volumes_index = np.array([el_handle.empty() for el_handle in check_volumes]).astype(bool)
            volumes_on_skin_handles = self.range_index(np.bitwise_not(external_volumes_index),self.all_volumes)
        elif self.dimension == 2:
            edges_on_skin_handles = self.bridge_adjacencies(self.all_edges,self.dimension)
            nodes_on_skin_handles = self.access_handle(edges_on_skin_handles)
            nodes_in_faces = ([self.mb.get_adjacencies(el_handle,0) for el_handle in self.all_faces])
            check_faces= ([rng.intersect(el_handle,nodes_on_skin_handles) for el_handle in nodes_in_faces])
            external_faces_index = np.array([el_handle.empty() for el_handle in check_faces]).astype(bool)
            faces_on_skin_handles = self.range_index(np.bitwise_not(external_faces_index),self.all_faces)
            volumes_on_skin_handles = rng.Range()
        return [nodes_on_skin_handles, edges_on_skin_handles, faces_on_skin_handles, volumes_on_skin_handles]
    def skinner_operation(self):
        self.skin = sk.Skinner(self.mb)
        if self.dimension == 3:
            faces_on_skin_handles  = self.skin.find_skin(0, self.all_volumes)
            edges_on_skin_handles = self.mtu.get_bridge_adjacencies(faces_on_skin_handles, 2, 1)
            nodes_on_skin_handles = self.mtu.get_bridge_adjacencies(faces_on_skin_handles, 2, 0)
            volumes_on_skin_handles = rng.intersect(self.mtu.get_bridge_adjacencies(faces_on_skin_handles, 0, 3), self.all_volumes)
        elif self.dimension == 2:
            edges_on_skin_handles = self.skin.find_skin(0, self.all_faces)
            nodes_on_skin_handles = self.mtu.get_bridge_adjacencies(edges_on_skin_handles, 1, 0)
            faces_on_skin_handles = rng.intersect(self.mtu.get_bridge_adjacencies(edges_on_skin_handles, 0, 2), self.all_faces)
            volumes_on_skin_handles = self.mtu.get_bridge_adjacencies(edges_on_skin_handles, 0, 3) #empty
        return [nodes_on_skin_handles, edges_on_skin_handles, faces_on_skin_handles, volumes_on_skin_handles]

    def bridge_adjacencies(self, handle, dim):
        # lacks support for indexing with multiple numbers
        if dim == 3:
            all_bridge = [self.mtu.get_bridge_adjacencies(el, 2, 3) for el in handle]
        else:
            all_bridge = [self.mtu.get_bridge_adjacencies(el, 1, 2) for el in handle]
        inside_meshset = self.mb.get_entities_by_handle(self.root_set)
        all_brige_in_meshset = np.array([rng.intersect(el_handle, inside_meshset) for el_handle in all_bridge])
        size_brige = np.array([len(el_handle) for el_handle in all_brige_in_meshset])
        handles = np.asarray(handle)[size_brige == 1].astype("uint")
        return rng.Range(handles)
