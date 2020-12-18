"""
Module for implementation of multiscale mesh and CoarseVolumes objects functionalities
"""
from . finescaleMesh import FineScaleMesh
from ..msCoarseningLib import algoritmo
from ..msCoarseningLib.partitionTools import partitionManager
from . serialization import IMPRESSPickler, IMPRESSUnpickler
from . meshComponents import MoabVariable
from . mscorePymoab import MsCoreMoab
from . meshComponentsMS import MoabVariableMS,  MeshEntitiesMS
from ..meshHandle.configTools.configClass import variableInit, coarseningInit
from pymoab import core, types, rng
import numpy as np
import pickle
from scipy.sparse import lil_matrix

print('Initializing Finescale Mesh for Multiscale Methods')

class FineScaleMeshMS(FineScaleMesh):
    def __init__(self, mesh_file, dim=3, var_config=None, load=False):
        self.var_config = var_config
        super().__init__(mesh_file, dim, load=load)
        print("Creating Coarse Grid")
        self.coarse = MultiscaleCoarseGrid(self, var_config, load = load)
        self.enhance_entities()

    def __getitem__ (self, key):
        if not isinstance(key, int):
            raise ValueError("Invalid key type provided")
        return self.coarse.elements[key]

    def enhance_entities(self):
        for i, el in zip(range(len(self.coarse.elements)), self.coarse.elements):
            el(i, self.coarse)

    def init_entities(self):
        self.nodes = MeshEntitiesMS(self.core, entity_type= "node")
        self.edges = MeshEntitiesMS(self.core, entity_type= "edges")
        self.faces = MeshEntitiesMS(self.core, entity_type= "faces")
        if self.dim == 3:
            self.volumes = MeshEntitiesMS(self.core, entity_type="volumes")

    def save_variables(self, name_file):
        self.core.mb.write_file('saves/' + name_file + '.h5m')
        with open('saves/' + name_file + '.imp', 'wb') as fp:
            pickle.dump([(tags.name_tag, tags.var_type, tags.data_size, tags.data_format, tags.data_density) for tags in self.var_handle_list], fp)
        for elements in self.coarse.elements:
            elements.save_variables(name_file)

    def load_variables(self):
        self.var_handle_list = []
        with open(self.mesh_file.split('.')[0]+'.imp', 'rb') as fp:
            tag_list = pickle.load(fp)
        for tags in tag_list:
            self.create_variable(name_tag = tags[0], var_type = tags[1], data_size = tags[2], data_format = tags[3], data_density = tags[4], create = False)

    def create_variable( self, name_tag, var_type= "volumes", data_size=1, data_format= "float", data_density= "sparse",
                 entity_index=None,  level=0, coarse_num=0, create = True):
        var = MoabVariableMS(self.core, data_size = data_size, var_type = var_type, data_format = data_format, name_tag = name_tag, data_density = data_density, entity_index = entity_index, level = level, coarse_num = coarse_num, create = create)
        exec(f'self.{name_tag} = var')
        self.var_handle_list.append(var)
        return var

    def to_moab(self):
        for vars in self.var_handle_list:
            vars.to_moab()
        for elements in self.coarse.elements:
            elements.to_moab()

    def to_numpy(self):
        for vars in self.var_handle_list:
            vars.to_numpy()
        for elements in self.coarse.elements:
            elements.to_numpy()

    def init_variables(self):
        self.var_handle_list = []
        if self.var_config is None:
            self.var_config = variableInit()
        for command in self.var_config.get_var(self.core.level):
            exec(command)

    def init_partition(self):
        coarse_config = coarseningInit()

        partitioner = partitionManager(self, coarse_config)
        [partition_tag, _] = partitioner()

        if isinstance(partition_tag, str) and partition_tag == 'parallel':
            return self.init_partition_parallel()
        else:
            if self.dim == 2:
                partition_moab = MoabVariable(self.core, data_size=1,
                                        var_type="faces", data_format="int",
                                        name_tag="Partition", data_density="dense")
            elif self.dim == 3:
                partition_moab = MoabVariable(self.core, data_size=1,
                                        var_type="volumes", data_format="int",
                                        name_tag="Partition", data_density="dense")
            partition_moab[:] = partition_tag
            return partition_moab

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
    def __init__(self, father_core, dim, i, coarse_vec, var_config=None , load = False, mesh_file = None):
        self.var_config = var_config
        self.dim = dim
        self.mesh_file = mesh_file
        self.level = father_core.level + 1
        self.coarse_num = i
        print("Level {0} - Volume {1}".format(self.level,self.coarse_num))
        self.core = MsCoreMoab(father_core, i, coarse_vec)
        self.init_entities()
        if not load:
            self.init_variables()
            self.init_coarse_variables()
        else:
            self.load_variables()
        self.macro_dim()

    def init_variables(self):
        self.var_handle_list = []
        if self.var_config is None:
            self.var_config = variableInit()
        for command in self.var_config.get_var(self.core.level, self.coarse_num):
            exec(command)

    def save_variables(self, name_file):
        name = self.core.id_name
        name = name[(name.find("ID") + 3):]
        with open('saves/'+name_file+name+'.imp', 'wb') as fp:
            pickle.dump([(tags.name_tag, tags.var_type, tags.data_size, tags.data_format, tags.data_density) for tags in self.var_handle_list], fp)

    def load_variables(self):
        self.var_handle_list = []
        name = self.core.id_name
        name = name[(name.find("ID") + 3):]
        with open(self.mesh_file.split('.')[0]+name+'.imp', 'rb') as fp:
            tag_list = pickle.load(fp)
        for tags in tag_list:
            self.create_variable(name_tag = tags[0], var_type = tags[1], data_size = tags[2], data_format = tags[3], data_density = tags[4], level = self.level, coarse_num = self.coarse_num, create = False, suffix = name)
        return

    def create_variable(self, name_tag, var_type="volumes", data_size=1, data_format="float", data_density="sparse",
                 entity_index=None,  level=0, coarse_num=0, create = True, suffix = None):
        var = MoabVariableMS(self.core, data_size = data_size, var_type = var_type, data_format = data_format, name_tag = name_tag, data_density = data_density, entity_index = entity_index, level = level, coarse_num = coarse_num, create = create)
        if suffix is not None:
            name_tag = name_tag.replace(suffix, '')
        exec(f'self.{name_tag} = var')
        self.var_handle_list.append(var)
        return var

    def to_moab(self):
        for variables in self.var_handle_list:
            variables.to_moab()

    def to_numpy(self):
        for variables in self.var_handle_list:
            variables.to_numpy()

    def __call__(self, i, general):
        self.nodes.enhance(i, general)
        self.edges.enhance(i, general)
        self.faces.enhance(i, general)
        if self.dim == 3:
            self.volumes.enhance(i, general)

    def init_coarse_variables(self):
        pass

class GetCoarseItem(object):
    def __init__(self, adj, tag, dic):
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
            return self.fun(self.tag, self.dic[item].get_array(), flat=True).astype(np.int64)
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
    def __init__(self, M, var_config, load = False):
        self.mb = M.core.mb
        self.M = M
        self.partition = M.init_partition()
        self.elements = [CoarseVolume(M.core, M.dim, i,
                        self.partition[:].ravel() == i, var_config)
                        for i in range(self.partition[:].max()+1 )]
        self.num_coarse = len(self.elements)
        self.num = {"nodes": 0, "node": 0, "edges": 1, "edge": 1, "faces": 2,
                    "face": 2, "volumes": 3, "volume": 3,0: 0, 1: 1, 2: 2, 3: 3}

        self.local_volumes_tag = [volume.core.handleDic[volume.core.id_name] for volume in self.elements]
        self.father_tag = M.core.handleDic[M.core.id_name]
        self.global_tag = M.core.handleDic["GLOBAL_ID"]
        self._all_volumes = M.core.all_volumes
        self._all_faces = M.core.all_faces
        self._all_edges = M.core.all_edges
        self._all_nodes = M.core.all_nodes
        if self.M.dim == 2:
            self.find_coarse_neighbours_2d()
            self.interfaces_edges = GetCoarseItem(self.mb.tag_get_data, self.father_tag, self._edges)
            self.interfaces_nodes = GetCoarseItem(self.mb.tag_get_data, self.father_tag, self._nodes)
            self.iface_coarse_neighbors = self._internal_edges(M)
        elif self.M.dim == 3:
            self.find_coarse_neighbours_3d()
            self.interfaces_faces = GetCoarseItem(self.mb.tag_get_data, self.father_tag, self._faces)
            self.interfaces_edges = GetCoarseItem(self.mb.tag_get_data, self.father_tag, self._edges)
            self.interfaces_nodes = GetCoarseItem(self.mb.tag_get_data, self.father_tag, self._nodes)
            self.iface_coarse_neighbors = self._internal_faces(M)

    def _internal_faces(self, M):
        faces = np.array([self._faces[el][0] for el in range (len(self._faces))], dtype=np.uint64)
        faces = self.mb.tag_get_data(self.father_tag, faces, flat=True)
        partition = self.partition[:].ravel()
        external = faces[self.num_internal_faces:]
        external_volumes = M.faces.bridge_adjacencies(external, interface="faces",target="volumes")
        ext_neigh = np.zeros((external_volumes.shape[0],2))
        ext_neigh[:,0], ext_neigh[:,1] = partition[external_volumes].ravel(), partition[external_volumes].ravel()
        if self.num_internal_faces == 0:
            return ext_neigh
        internal = faces[0:self.num_internal_faces]
        internal_volumes = M.faces.bridge_adjacencies(internal, interface="faces",target="volumes")
        int_neigh = np.vstack((partition[internal_volumes[:,0]],partition[internal_volumes[:,1]])).T
        return np.vstack((int_neigh,ext_neigh)).astype(np.int64)
    
    def _internal_edges(self, M):
        edges = np.array([self._edges[el][0] for el in range (len(self._edges))], dtype=np.uint64)
        edges = self.mb.tag_get_data(self.father_tag, edges, flat=True)
        partition = self.partition[:].ravel()
        external = edges[self.num_internal_edges:]
        external_faces = M.edges.bridge_adjacencies(external, interface="edges",target="faces")
        ext_neigh = np.zeros((external_faces.shape[0],2))
        ext_neigh[:,0], ext_neigh[:,1] = partition[external_faces].ravel(), partition[external_faces].ravel()
        if self.num_internal_edges == 0:
            return ext_neigh
        internal = edges[0:self.num_internal_edges]
        internal_faces = M.faces.bridge_adjacencies(internal, interface="edges",target="faces")
        int_neigh = np.vstack((partition[internal_faces[:,0]],partition[internal_faces[:,1]])).T
        return np.vstack((int_neigh,ext_neigh)).astype(np.int64)

    def find_coarse_neighbours_2d(self):
        self.all_nodes_neighbors = rng.Range()
        self.all_edges_neighbors = rng.Range()
        self.all_faces_neighbors = rng.Range()

        self._faces_neighbors = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.uint32)
        self._edges_neighbors = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.uint32)
        self._nodes_neighbors = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.uint32)

        self.faces_connectivities = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.bool)
        self.edges_connectivities = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.bool)
        self.nodes_connectivities = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.bool)

        edges_array = self.M.core.internal_edges.get_array()
        adj_array = self.mb.get_ord_adjacencies(edges_array, 2)[0]

        tg = self.mb.tag_get_handle('Partition')
        parts = self.mb.tag_get_data(tg, adj_array).reshape(-1,2)
        
        boundaries = self.M.core.boundary_edges.get_array()
        boundary_face = self.M.core.mb.get_ord_adjacencies(boundaries, 2)[0]
        self.all_faces_neighbors.insert(boundary_face)
        self.all_edges_neighbors.insert(boundaries)
        boundary_parts = self.mb.tag_get_data(tg, boundary_face, flat=True)

        indx = np.where(parts[:,0] != parts[:,1])[0]
        parts = parts[indx]
        inters_edges = edges_array[indx]
        self._edges, self.num_internal_edges = self.M.core.mb.get_interface_faces2(
            self.edges_connectivities, parts, inters_edges, boundaries, boundary_parts, 
            self.num_coarse, self._edges_neighbors)
        
        print('Matrix of coarse edges adjacencies created')

        if inters_edges.size == 0:
            inters_nodes = np.array([], dtype=np.uint64)
            indx = np.array([], dtype=np.int64)
            coarse_jagged = np.array([], dtype=np.uint64)
        else:
            self.all_edges_neighbors.insert(inters_edges)
            self.all_faces_neighbors.insert(adj_array.reshape(-1,2)[indx].ravel())
            inters_nodes = np.unique(self.mb.get_ord_adjacencies(inters_edges, 0)[0])
            self.all_nodes_neighbors.insert(inters_nodes)
            aux_tuple = self.M.core.mb.get_ord_adjacencies(inters_nodes, 2)
            temp_jagged = np.delete(np.array(np.split(aux_tuple[0], aux_tuple[1]), dtype=object), -1)
            jagged_index = np.array([temp_jagged[i].size 
                for i in range(temp_jagged.shape[0])], dtype=np.int32)
            jagged_index = np.cumsum(jagged_index, dtype=np.int32)[:-1]
            coarse_array = self.M.core.mb.tag_get_data(tg, np.concatenate(temp_jagged), flat=True)
            coarse_jagged = np.array(np.split(coarse_array, jagged_index), dtype=object)
            coarse_jagged = np.array([np.unique(coarse_jagged[i]) 
                for i in range(coarse_jagged.shape[0])], dtype=object)
            indx = np.array([coarse_jagged[i].size > 2 for i in range(coarse_jagged.shape[0])])
        
        boundaries = self.M.core.boundary_nodes.get_array()
        self.all_nodes_neighbors.insert(boundaries)
        aux_tuple = self.M.core.mb.get_ord_adjacencies(boundaries, 2)
        temp_jagged = np.delete(np.array(np.split(aux_tuple[0], aux_tuple[1]), dtype=object), -1)
        jagged_index = np.array([temp_jagged[i].size for i in range(temp_jagged.shape[0])], dtype=np.int32)
        jagged_index = np.cumsum(jagged_index, dtype=np.int32)[:-1]

        boundary_parts = self.M.core.mb.tag_get_data(tg, np.concatenate(temp_jagged), flat=True)
        boundary_parts = np.array(np.split(boundary_parts, jagged_index), dtype=object)
        boundary_parts = np.array([np.unique(boundary_parts[i]).astype(np.int32) \
                                for i in range(boundary_parts.shape[0])], dtype=np.object)
        
        self._nodes, self.num_internal_nodes = self.M.core.mb.get_interface_entities2(
            self.nodes_connectivities, inters_nodes, coarse_jagged, indx, 
            boundaries, boundary_parts, self.num_coarse, self._nodes_neighbors)
        
        print('Matrix of coarse nodes adjacencies created')
        
        self._faces_neighbors = self._faces_neighbors.tocsr()
        self._edges_neighbors = self._edges_neighbors.tocsr()
        self._nodes_neighbors = self._nodes_neighbors.tocsr()
        self.faces_connectivities = self.faces_connectivities.tocsr()
        self.edges_connectivities = self.edges_connectivities.tocsr()
        self.nodes_connectivities = self.nodes_connectivities.tocsr()
        self.connectivities = (self.nodes_connectivities, self.edges_connectivities, self.faces_connectivities)
    
    def find_coarse_neighbours_3d(self):
        self.all_nodes_neighbors = rng.Range()
        self.all_edges_neighbors = rng.Range()
        self.all_faces_neighbors = rng.Range()
        self.all_volumes_neighbors = rng.Range()

        self._faces_neighbors  = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.uint32)
        self._edges_neighbors  = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.uint32)
        self._nodes_neighbors  = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.uint32)

        self.faces_connectivities = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.bool)
        self.edges_connectivities = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.bool)
        self.nodes_connectivities = lil_matrix((self.num_coarse,self.num_coarse+1), dtype=np.bool)

        faces_array = self.M.core.internal_faces.get_array()
        adj_array = self.mb.get_ord_adjacencies(faces_array, 3)[0]

        tg = self.mb.tag_get_handle('Partition')
        parts = self.mb.tag_get_data(tg, adj_array).reshape(-1,2)

        boundaries = self.M.core.boundary_faces.get_array()
        boundary_vol = self.M.core.mb.get_ord_adjacencies(boundaries, 3)[0]
        self.all_volumes_neighbors.insert(boundary_vol)
        self.all_faces_neighbors.insert(boundaries)
        boundary_parts = self.mb.tag_get_data(tg, boundary_vol, flat=True)

        indx = np.where(parts[:,0]!=parts[:,1])[0]
        parts = parts[indx]
        inters_faces = faces_array[indx]
        self._faces, self.num_internal_faces = self.M.core.mb.get_interface_faces2(
            self.faces_connectivities, parts, inters_faces, boundaries, boundary_parts, 
            self.num_coarse, self._faces_neighbors)
        
        print('Matrix of coarse faces adjacencies created')

        if inters_faces.size == 0:
            inters_edges = np.array([], dtype=np.uint64)
            indx = np.array([], dtype=np.int64)
            coarse_jagged = np.array([], dtype=np.uint64)
        else:
            self.all_faces_neighbors.insert(inters_faces)
            self.all_volumes_neighbors.insert(adj_array.reshape(-1,2)[indx].ravel())
            inters_edges = np.unique(self.mb.get_ord_adjacencies(inters_faces, 1)[0])
            self.all_edges_neighbors.insert(inters_edges)
            aux_tuple = self.M.core.mb.get_ord_adjacencies(inters_edges, 3)
            temp_jagged = np.delete(np.array(np.split(aux_tuple[0], aux_tuple[1]), dtype=object), -1)
            jagged_index = np.array([temp_jagged[i].size 
                for i in range(temp_jagged.shape[0])], dtype=np.int32)
            jagged_index = np.cumsum(jagged_index, dtype=np.int32)[:-1]
            coarse_array = self.M.core.mb.tag_get_data(tg, np.concatenate(temp_jagged), flat=True)
            coarse_jagged = np.array(np.split(coarse_array, jagged_index), dtype=object)
            coarse_jagged = np.array([np.unique(coarse_jagged[i]) 
                for i in range(coarse_jagged.shape[0])], dtype=object)
            indx = np.array([coarse_jagged[i].size > 2 for i in range(coarse_jagged.shape[0])])

        boundaries = self.M.core.boundary_edges.get_array()
        self.all_edges_neighbors.insert(boundaries)
        aux_tuple = self.M.core.mb.get_ord_adjacencies(boundaries, 3)
        temp_jagged = np.delete(np.array(np.split(aux_tuple[0], aux_tuple[1]), dtype=object), -1)
        jagged_index = np.array([temp_jagged[i].size for i in range(temp_jagged.shape[0])], dtype=np.int32)
        jagged_index = np.cumsum(jagged_index, dtype=np.int32)[:-1]

        boundary_parts = self.M.core.mb.tag_get_data(tg, np.concatenate(temp_jagged), flat=True)
        boundary_parts = np.array(np.split(boundary_parts, jagged_index), dtype=object)
        boundary_parts = np.array([np.unique(boundary_parts[i]).astype(np.int32) \
                                for i in range(boundary_parts.shape[0])], dtype=np.object)
        
        self._edges, self.num_internal_edges = self.M.core.mb.get_interface_entities2(
            self.edges_connectivities, inters_edges, coarse_jagged, indx, 
            boundaries, boundary_parts, self.num_coarse, self._edges_neighbors)
        
        print('Matrix of coarse edges adjacencies created')

        if inters_faces.size == 0:
            inters_nodes = np.array([], dtype=np.uint64)
            indx = np.array([], dtype=np.int64)
            coarse_jagged = np.array([], dtype=np.uint64)
        else:
            inters_nodes = np.unique(self.mb.get_ord_adjacencies(inters_faces, 0)[0])
            self.all_nodes_neighbors.insert(inters_nodes)
            aux_tuple = self.M.core.mb.get_ord_adjacencies(inters_nodes, 3)
            temp_jagged = np.delete(np.array(np.split(aux_tuple[0], aux_tuple[1]), dtype=object), -1)
            jagged_index = np.array([temp_jagged[i].size for i in range(temp_jagged.shape[0])], dtype=np.int32)
            jagged_index = np.cumsum(jagged_index, dtype=np.int32)[:-1]
            coarse_array = self.M.core.mb.tag_get_data(tg, np.concatenate(temp_jagged), flat=True)
            coarse_jagged = np.array(np.split(coarse_array, jagged_index), dtype=object)
            coarse_jagged = np.array([np.unique(coarse_jagged[i]) 
                for i in range(coarse_jagged.shape[0])], dtype=object)
            indx = np.array([coarse_jagged[i].size>2 for i in range(coarse_jagged.shape[0])])

        boundaries = self.M.core.boundary_nodes.get_array()
        self.all_nodes_neighbors.insert(boundaries)
        aux_tuple = self.M.core.mb.get_ord_adjacencies(boundaries, 3)

        temp_jagged = np.delete(np.array(np.split(aux_tuple[0], aux_tuple[1]), dtype=object), -1)
        jagged_index = np.array([temp_jagged[i].size for i in range(temp_jagged.shape[0])], dtype=np.int32)
        jagged_index = np.cumsum(jagged_index, dtype=np.int32)[:-1]
        boundary_parts = self.M.core.mb.tag_get_data(tg, np.concatenate(temp_jagged), flat=True)
        boundary_parts = np.array(np.split(boundary_parts, jagged_index), dtype=object)
        boundary_parts = np.array([np.unique(boundary_parts[i]) 
            for i in range(boundary_parts.shape[0])], dtype=object)
        self._nodes, self.num_internal_nodes = self.M.core.mb.get_interface_entities2(
            self.nodes_connectivities, inters_nodes, coarse_jagged, indx, boundaries, 
            boundary_parts, self.num_coarse, self._nodes_neighbors)
        
        print('Matrix of coarse nodes adjacencies created')

        self._faces_neighbors = self._faces_neighbors.tocsr()
        self._edges_neighbors = self._edges_neighbors.tocsr()
        self._nodes_neighbors = self._nodes_neighbors.tocsr()
        self.faces_connectivities = self.faces_connectivities.tocsr()
        self.edges_connectivities = self.edges_connectivities.tocsr()
        self.nodes_connectivities = self.nodes_connectivities.tocsr()
        self.connectivities = (self.nodes_connectivities, self.edges_connectivities, self.faces_connectivities)

    def iface_neighbors(self, x):
        tmp = self._faces_neighbors[x].A[0]
        indx = np.nonzero(tmp)[0]
        return indx, (tmp[indx]-1).astype(np.int64)

    def iedge_neighbors(self, x):
        tmp = self._edges_neighbors[x].A[0]
        indx = np.nonzero(tmp)[0]
        return indx, (tmp[indx]-1).astype(np.int64)

    def inode_neighbors(self, x):
        tmp = self._nodes_neighbors[x].A[0]
        indx = np.nonzero(tmp)[0]
        return indx, (tmp[indx]-1).astype(np.int64)

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
              return self.mb.tag_get_data(self.father_tag, 
                self._nodes[self._nodes_neighbors[x,y]-1].get_array(), flat=True).astype(np.int64)
          elif flag == 1:
              return self.mb.tag_get_data(self.father_tag, 
                self._edges[self._edges_neighbors[x,y]-1].get_array(), flat=True).astype(np.int64)
          elif flag == 2:
              return self.mb.tag_get_data(self.father_tag, 
                self._faces[self._faces_neighbors[x,y]-1].get_array(), flat=True).astype(np.int64)

    @property
    def all_interface_nodes(self):
        return self.mb.tag_get_data(self.father_tag, 
            self.all_nodes_neighbors.get_array(), flat=True).astype(np.int64)

    @property
    def all_interface_edges(self):
        return self.mb.tag_get_data(self.father_tag, 
            self.all_edges_neighbors.get_array(), flat=True).astype(np.int64)

    @property
    def all_interface_faces(self):
        return self.mb.tag_get_data(self.father_tag, 
            self.all_faces_neighbors.get_array(), flat=True).astype(np.int64)

    @property
    def all_interface_volumes(self):
        return self.mb.tag_get_data(self.father_tag, 
            self.all_volumes_neighbors.get_array(), flat=True).astype(np.int64)

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
            range_el = self._all_volumes
        if index_vec.size > 0:
            range_el = self.range_index(index_vec,range_el)
        try:
            handle_tag = self.M.core.handleDic[name_tag]
            return self.mb.tag_get_data(handle_tag, range_el)
        except KeyError:
            print("Tag not found")

    def range_index(self, vec_index, range_handle = None):
        if range_handle is None:
            range_handle = self._all_volumes
        if vec_index.dtype == "bool":
            vec = np.where(vec_index)[0]
        else:
            vec = vec_index.astype("uint")
        handles = np.asarray(range_handle)[vec.astype("uint")].astype("uint")
        return rng.Range(handles)
