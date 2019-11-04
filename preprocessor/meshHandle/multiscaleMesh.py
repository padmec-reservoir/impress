"""
Module for implementation of multiscale mesh and CoarseVolumes objects functionalities
"""

import pdb
from . finescaleMesh import FineScaleMesh
from ..msCoarseningLib import algoritmo
from ..msCoarseningLib.partitionTools import partitionManager
from . meshComponents import MoabVariable
from . mscorePymoab import MsCoreMoab
from . meshComponentsMS import MoabVariableMS,  MeshEntitiesMS
from ..meshHandle.configTools.configClass import variableInit, coarseningInit
from pymoab import core, types, rng
import numpy as np



print('Initializing Finescale Mesh for Multiscale Methods')


class FineScaleMeshMS(FineScaleMesh):
    def __init__(self, mesh_file, dim=3, var_config=None):
        self.var_config = var_config
        super().__init__(mesh_file, dim)
        print("Creating Coarse Grid")
        self.coarse = MultiscaleCoarseGrid(self, var_config)
        self.enhance_entities()

    def __getitem__ (self, key):
        if not isinstance(key, int):
            raise ValueError("Invalid key type provided")
        return self.coarse.elements[key]

    def enhance_entities(self):
        for i,el in zip(range(len(self.coarse.elements)), self.coarse.elements):
            el(i,self.coarse)

    def init_entities(self):
        self.nodes = MeshEntitiesMS(self.core, entity_type="node")
        self.edges = MeshEntitiesMS(self.core, entity_type="edges")
        self.faces = MeshEntitiesMS(self.core, entity_type="faces")
        if self.dim == 3:
            self.volumes = MeshEntitiesMS(self.core, entity_type="volumes")

    def init_variables(self):
        if self.var_config is None:
            self.var_config = variableInit()
        for command in self.var_config.get_var(self.core.level):
            exec(command)

    def init_partition(self):
        coarse_config = coarseningInit()
        partitioner = partitionManager(self, coarse_config)
        [partition, coarse_center] = partitioner()
        if isinstance(partition, str) and partition == 'parallel':
            return self.init_partition_parallel()
        else:
            return partition

    def init_partition_parallel(self):
        if self.dim == 3:
            partition = MoabVariable(self.core, data_size=1,
                                     var_type="volumes", data_format="int",
                                     name_tag="Parallel", data_density="sparse")
        elif self.dim == 2:
            partition = MoabVariable(self.core, data_size=1, var_type="faces",
                                     data_format="int", name_tag="Parallel",
                                     data_density="sparse")
        return partition


class CoarseVolume(FineScaleMeshMS):
    def __init__(self, father_core, dim, i, coarse_vec, var_config=None):
        self.var_config = var_config
        self.dim = dim
        self.level = father_core.level + 1
        self.coarse_num = i
        self.coarse_vec = coarse_vec
        print("Level {0} - Volume {1}".format(self.level,self.coarse_num))
        self.core = MsCoreMoab(father_core, i, coarse_vec)
        self.init_entities()
        self.init_variables()
        self.init_coarse_variables()
        self.macro_dim()

    def init_variables(self):
        if self.var_config is None:
            self.var_config = variableInit()
        for command in self.var_config.get_var(self.core.level, self.coarse_num):
            exec(command)

    def __call__(self, i, general):
        self.nodes.enhance(i, general)
        self.edges.enhance(i, general)
        self.faces.enhance(i, general)
        if self.dim == 3:
            self.volumes.enhance(i, general)

    def init_coarse_variables(self):
        pass


class GetCoarseItem(object):
    def __init__(self, adj,tag, dic):
        self.fun = adj
        self.tag = tag
        self.dic = dic

    def __len__(self):
        return len(self.dic)

    def __call__(self, item):
        tmp = self.dic[item]
        el_list = rng.Range()
        for e in tmp:
            el_list.insert(e[0])
        return self.fun(self.tag, el_list)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.fun(self.tag, self.dic[item]).ravel()
        elif isinstance(item, slice):
            start = item.start
            step = item.step
            stop = item.stop
            if step == None:
                step = 1
            if start == None:
                start = 0
            if stop == None:
                stop = len(self.dic)
            array = np.array(range(start, stop, step))
            s = np.array([])
            for el in array:
                s = np.concatenate((s, self.__getitem__(int(el))))
            return s


class MultiscaleCoarseGrid(object):
    def __init__(self, M, var_config):
        self.mb = M.core.mb
        self.partition = M.init_partition()
        self.elements = [CoarseVolume(M.core, M.dim, i, self.partition[:].ravel() == i, var_config) for i in range(self.partition[:].max()+1 )]
        self.num_coarse = len(self.elements)
        self.num = {"nodes": 0, "node": 0, "edges": 1, "edge": 1, "faces": 2, "face": 2, "volumes": 3, "volume": 3,
                             0: 0, 1: 1, 2: 2, 3: 3}
        self.local_volumes_tag = [volume.core.handleDic[volume.core.id_name] for volume in self.elements]
        self.father_tag = M.core.handleDic[M.core.id_name]
        self.global_tag = M.core.handleDic["GLOBAL_ID"]
        self._all_volumes = M.core.all_volumes
        self._all_faces = M.core.all_faces
        self._all_edges = M.core.all_edges
        self._all_nodes = M.core.all_nodes
        self.find_coarse_neighbours()
        self.interfaces_faces = GetCoarseItem(self.mb.tag_get_data, self.father_tag, self._faces)
        self.interfaces_edges = GetCoarseItem(self.mb.tag_get_data, self.father_tag, self._edges)
        self.interfaces_nodes = GetCoarseItem(self.mb.tag_get_data, self.father_tag, self._nodes)
        self.iface_coarse_neighbors = self._internal_faces(M)

    def _internal_faces(self, M):
        #faces = np.array([self.mb.tag_get_data(self.father_tag,el[0]).ravel() for el in self._faces]).ravel()
        #faces = [self.interfaces_faces[int(el)] for el in range(len(self.interfaces_faces))]
        faces = []
        for el in range(len(self.interfaces_faces)):
            faces.append(self.interfaces_faces[el][0])
        partition = self.partition[:].ravel()
        internal = faces[0:self.num_internal_faces]
        external = faces[self.num_internal_faces:]
        internal_volumes = M.faces.bridge_adjacencies(internal, interface="faces",target="volumes")
        external_volumes = M.faces.bridge_adjacencies(external, interface="faces",target="volumes")
        int_neigh = np.vstack((partition[internal_volumes[:,0]],partition[internal_volumes[:,1]])).T
        ext_neigh = np.zeros((external_volumes.shape[0],2))
        ext_neigh[:,0], ext_neigh[:,1] = partition[external_volumes].ravel(), partition[external_volumes].ravel()
        return np.vstack((int_neigh,ext_neigh)).astype("uint64")

    def find_coarse_neighbours(self):
        self.connectivities = np.zeros((self.num_coarse,self.num_coarse+1 ,3)).astype('bool')
        self._nodes_neighbors  = np.zeros((self.num_coarse,self.num_coarse+1), dtype = object)
        self._edges_neighbors  = np.zeros((self.num_coarse,self.num_coarse+1), dtype = object)
        self._faces_neighbors  = np.zeros((self.num_coarse,self.num_coarse+1), dtype = object)
        self._nodes_neighbors[:], self._edges_neighbors[:], self._faces_neighbors[:] = None , None, None
        self._nodes = list()
        self._faces = list()
        self._edges = list()
        self.all_nodes_neighbors = rng.Range()
        self.all_edges_neighbors = rng.Range()
        self.all_faces_neighbors = rng.Range()
        self.all_volumes_neighbors = rng.Range()
        #print(self.all_nodes_neighbors, self.all_edges_neighbors, self.all_faces_neighbors)
        # import pdb; pdb.set_trace()
        node_count, edge_count, face_count = 0, 0, 0
        for x in range(self.num_coarse):
            for y in range(x+1,self.num_coarse):
                node_intersect = rng.intersect(self.elements[x].core.boundary_nodes, self.elements[y].core.boundary_nodes)
                if not node_intersect.empty():
                    self._nodes.append(node_intersect)
                    #self._nodes = np.append(self._nodes,node_intersect)
                    self._nodes_neighbors[x,y], self._nodes_neighbors[y,x],= node_count ,node_count
                    self.connectivities[x, y, 0],self.connectivities[y, x, 0] = True, True
                    node_count += 1
                    [self.all_nodes_neighbors.insert(e) for e in node_intersect]
                edges_intersect = rng.intersect(self.elements[x].core.boundary_edges, self.elements[y].core.boundary_edges)
                if not edges_intersect.empty():
                    self._edges.append(edges_intersect)
                    # self._edges = np.append(self._edges,edges_intersect)
                    self._edges_neighbors[x,y], self._edges_neighbors[y,x]= edge_count ,edge_count
                    self.connectivities[x, y, 1], self.connectivities[y, x, 1] =  True, True
                    edge_count += 1
                    [self.all_edges_neighbors.insert(e) for e in edges_intersect]
                faces_intersect = rng.intersect(self.elements[x].core.boundary_faces, self.elements[y].core.boundary_faces)
                if not faces_intersect.empty():
                    self._faces.append(faces_intersect)
                    #self._faces = np.append(self._faces,faces_intersect)
                    self._faces_neighbors[x,y], self._faces_neighbors[y,x]= face_count ,face_count
                    self.connectivities[x, y, 2],self.connectivities[y, x, 2]  = True, True
                    face_count += 1
                    [self.all_faces_neighbors.insert(e) for e in faces_intersect]
        self.num_internal_nodes = node_count
        self.num_internal_edges = edge_count
        self.num_internal_faces = face_count

        for x in range(self.num_coarse):
            #  fix the interesection - second variable poorly choosen
            node_intersect = rng.subtract(self.elements[x].core.boundary_nodes, self.all_nodes_neighbors)
            if not node_intersect.empty():
                self._nodes.append(node_intersect)
                self._nodes_neighbors[x, -1] = node_count
                self.connectivities[x, -1, 0] = True
                node_count += 1
            edge_intersect = rng.subtract(self.elements[x].core.boundary_edges, self.all_edges_neighbors)
            if not edge_intersect.empty():
                self._edges.append(edge_intersect)
                self._edges_neighbors[x, -1] = edge_count
                self.connectivities[x, -1, 1] = True
                edge_count += 1
            face_intersect = rng.subtract(self.elements[x].core.boundary_faces, self.all_faces_neighbors)
            if not face_intersect.empty():
                self._faces.append(face_intersect)
                self._faces_neighbors[x, -1] = face_count
                self.connectivities[x, -1, 2] = True
                face_count += 1

    def iface_neighbors(self, x):
        tmp = -1* np.ones(self._faces_neighbors[x].shape)
        tag = self._faces_neighbors[x] != None
        tmp[tag] = self._faces_neighbors[x,tag]
        #import pdb; pdb.set_trace()
        indices = np.where(tmp >= 0)[0]
        return indices, tmp[indices].astype(int)

    def iedge_neighbors(self, x):
        tmp = -1* np.ones(self._edges_neighbors[x].shape)
        tag = self._edges_neighbors[x] != None
        tmp[tag] = self._edges_neighbors[x,tag]
        #import pdb; pdb.set_trace()
        indices = np.where(tmp >= 0)[0]
        return indices, tmp[indices].astype(int)

    def inode_neighbors(self, x):
        tmp = -1* np.ones(self._nodes_neighbors[x].shape)
        tag = self._nodes_neighbors[x] != None
        tmp[tag] = self._nodes_neighbors[x,tag]
        #import pdb; pdb.set_trace()
        indices = np.where(tmp >= 0)[0]
        return indices, tmp[indices].astype(int)

    def father_to_local_id(self, vec_range,  element, target):
        flag = self.num[element]
        vec = self.create_range_vec(vec_range)
        if flag == 0:
            handle = self.range_index(vec, self._all_nodes)
        elif flag == 1:
            handle = self.range_index(vec, self._all_edges)
        elif flag == 2:
            handle = self.range_index(vec, self._all_faces)
        elif flag == 3:
            handle = self.range_index(vec, self._all_volumes)
        return self.mb.tag_get_data(self.local_volumes_tag[target],handle)

    def neighbours(self, x,y, element):
          flag = self.num[element]
          if flag == 0:
              return self.mb.tag_get_data(self.father_tag, self._nodes[self.nodes_neighbors[x,y]])
          elif flag == 1:
              return self.mb.tag_get_data(self.father_tag, self._edges[self.edges_neighbors[x,y]])
          elif flag == 2:
              return self.mb.tag_get_data(self.father_tag, self._faces[self.faces_neighbors[x,y]])

    @property
    def all_interface_nodes(self):
        return self.mb.tag_get_data(self.father_tag, self.all_nodes_neighbors)

    @property
    def all_interface_edges(self):
        return self.mb.tag_get_data(self.father_tag, self.all_edges_neighbors)

    @property
    def all_interface_faces(self):
        return self.mb.tag_get_data(self.father_tag, self.all_faces_neighbors)

    @property
    def all_neighbors_volumes(self):
        return self.mb.tag_get_data(self.father_tag, self.all_volumes_neighbors)

    def create_range_vec(self, index):
        range_vec = None
        if isinstance(index, int) or isinstance(index, np.integer):
            range_vec = np.array([index]).astype("uint")
        elif isinstance(index, np.ndarray):
            if index.dtype == "bool":
                range_vec = np.where(index)[0]
            else:
                range_vec = index
        elif isinstance(index, slice):
            start = index.start
            stop = index.stop
            step = index.step
            if start is None:
                start = 0
            if stop is None:
                stop = len(self)
            if step is None:
                step = 1
            if start < 0:
                start = len(self) + start + 1
            if stop < 0:
                stop = len(self) + stop + 1
            range_vec = np.arange(start, stop, step).astype('uint')
        elif isinstance(index, list):
            range_vec = np.array(index)
        return range_vec

    def read_data(self, name_tag, index_vec = np.array([]), range_el = None):
        if range_el is None:
            range_el = self.all_volumes
        if index_vec.size > 0:
            range_el = self.range_index(index_vec,range_el)
        try:
            handle_tag = self.handleDic[name_tag]
            return self.mb.tag_get_data(handle_tag, range_el)
        except KeyError:
            print("Tag not found")

    def range_index(self, vec_index, range_handle = None):
        if range_handle is None:
            range_handle = self.all_volumes
        if vec_index.dtype == "bool":
            vec = np.where(vec_index)[0]
        else:
            vec = vec_index.astype("uint")
        handles = np.asarray(range_handle)[vec.astype("uint")].astype("uint")
        return rng.Range(handles)
