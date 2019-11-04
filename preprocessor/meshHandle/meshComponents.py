"""
Generator of mesh entities and tags
"""

import pdb
import numpy as np
from pymoab import types, rng, topo_util
from ..geoUtil import geoTools as gtool
#from preprocessor.meshHandle.configTools.configClasses import variableInit


class GetItem(object):
    def __init__(self, adj):
        self.fun = adj

    def __call__(self, item):
        return self.fun(item)

    def __getitem__(self, item):
        return self.fun(item)


class MeshEntities(object):
    def __init__(self, core, entity_type):
        self.mb = core.mb
        self.mtu = core.mtu
        self.meshset = core.root_set
        self.nodes = core.all_nodes
        self.num = {"nodes": 0, "node": 0, "edges": 1, "edge": 1, "faces": 2, "face": 2, "volumes": 3, "volume": 3,
                    0: 0, 1: 1, 2: 2, 3: 3}
        string = {0: "nodes", 1: "edges", 2: "faces", 3: "volumes"}
        entity_num = self.num[entity_type]
        list_type = [(core.all_nodes, core.internal_nodes, core.boundary_nodes),
                        (core.all_edges, core.internal_edges, core.boundary_edges),
                        (core.all_faces, core.internal_faces, core.boundary_faces),
                        (core.all_volumes, core.internal_volumes, core.boundary_volumes)]
        self.list_all = [list_type[0][0], list_type[1][0], list_type[2][0], list_type[3][0]]
        self.level = core.level
        if self.level == 0:
            self.id_name = "GLOBAL_ID"
            self.father_id_name = "GLOBAL_ID"
            self.id_name = "GLOBAL_ID"
        elif self.level == 1:
            self.father_id_name = core.father_core.id_name
            self.id_name = "LOCAL_ID_L" + str(core.level) + "-" + str(core.coarse_num)
        else:
            self.father_id_name = core.father_core.id_name
            self.id_name = self.father_id_name + str("L") + str(self.level) + "-" + str(self.coarse_num)

        (self.elements_handle, self.internal_elements, self.boundary_elements), self.vID = list_type[entity_num], entity_num
        # self.internal_elements_array = self.internal_elements.get_array()
        self.entity_type = string[entity_num]
        self.tag_handle = core.handleDic[self.id_name]
        self.global_handle = core.handleDic['GLOBAL_ID']
        self.father_handle = core.handleDic[self.father_id_name]
        if self.vID == 0:
            self.adjacencies = GetItem(self._adjacencies_for_nodes)
            self.coords =  GetItem(self._coords)
        else:
            self.adjacencies = GetItem(self._adjacencies)
            self.connectivities = GetItem(self._connectivities)
        self.classify_element = GetItem(self._classify_element)
        self.center = GetItem(self._center)
        # self.global_id = GetItem(self._global_id)
        # self.father_id = GetItem(self._father_id)
        if (self.vID == 1) & (core.dimension == 2):
            self.normal = GetItem(self._normal)
        elif (self.vID == 2) & (core.dimension == 3):
            self.normal = GetItem(self._normal)
        # initialize specific flag dic in accordance with type of the object create
        self.flag = {key: self.read(value[self.vID]) for key, value in core.flag_dic.items()
                     if value[self.vID].empty() is not True}
        # print("Mesh Entity type {0} successfully initialized".format(entity_type))

    def bridge_adjacencies(self, index, interface, target):
        if isinstance(index, np.int64) or isinstance(index, int):
            el_handle = np.array([self.elements_handle[index]])
        elif not isinstance(index, np.ndarray) and index is not None:
            el_handle = self.elements_handle[index].get_array()
        else:
            el_handle = self.elements_handle.get_array(index)
        if self.level==0:
            return self.mtu.get_ord_bridge_adjacencies(el_handle, self.num[interface], self.num[target], self.mb, self.tag_handle)
        return self.mtu.get_ord_bridge_adjacencies(el_handle, self.num[interface], self.num[target], self.mb, self.tag_handle, self.list_all[self.num[target]], self.level)

    def _coords(self, index):
        if isinstance(index, np.int64) or isinstance(index, int):
            el_handle = np.array([self.elements_handle[index]])
        elif not isinstance(index, np.ndarray) and index is not None:
            el_handle = self.nodes[index]
        else:
            el_handle = self.nodes.get_array(index)
        return np.reshape(self.mb.get_coords(el_handle),(-1,3))

    # def _global_id(self, index):
    #     range_vec = self.create_range_vec(index)
    #     elements_handle = self.range_index(range_vec)
    #     return self.mb.tag_get_data(self.global_handle, elements_handle).ravel()
    # def _father_id(self, index):
    #     range_vec = self.create_range_vec(index)
    #     elements_handle = self.range_index(range_vec)
    #     return self.mb.tag_get_data(self.father_handle, elements_handle).ravel()

    def _adjacencies_for_nodes(self, index):
        return index

    def _adjacencies(self, index, flag_nodes=False, dim_tag=None):
        if dim_tag is None:
            if not flag_nodes:
                dim_tag = self.vID - 1
            else:
                dim_tag = 0
        if isinstance(index, np.int64) or isinstance(index, int):
            el_handle = np.array([self.elements_handle[index]])
        elif not isinstance(index, np.ndarray) and index is not None:
            el_handle = self.elements_handle[index].get_array()
        else:
            el_handle = self.elements_handle.get_array(index)
        return self.mb.get_ord_adjacencies(el_handle, dim_tag, tag_handle = self.tag_handle)

    def _center(self,index):
        if self.vID == 0:
            centers = self._coords(index)
            return centers
        elif self.vID == 1:
            edges_adj = self.connectivities[index]
            if edges_adj.ndim == 1:
                return 0.5*(self._coords(edges_adj[0]) + self._coords(edges_adj[1]))
            v0 = np.array([edges_adj[i][0] for i in range (edges_adj.shape[0])])
            v1 = np.array([edges_adj[i][1] for i in range (edges_adj.shape[0])])
            centers = 0.5 * (self._coords(v0) + self._coords(v1))
            return centers
        elif self.vID == 2 or self.vID == 3:
            if isinstance(index, np.int64) or isinstance(index, int):
                el_handle = np.array([self.elements_handle[index]])
            elif not isinstance(index, np.ndarray) and index is not None:
                el_handle = self.elements_handle[index].get_array()
            else:
                el_handle = self.elements_handle.get_array(index)
            adj = self.mb.get_ord_connectivity(el_handle, tag_opt = False)
            if adj.ndim==1:
                return gtool.get_average(np.reshape(self.mb.get_coords(adj),(-1,3)))
            centers = np.empty((adj.shape[0],3))
            for i in range(adj.shape[0]):
                centers[i] = gtool.get_average(np.reshape(self.mb.get_coords(adj[i]),(-1,3)))
                #pdb.set_trace()
            return centers
        return None

    def _normal(self,index):

        #normal_vec = np.zeros(( np.shape(range_vec)[0],3 ))
        adj = self.connectivities[index]
        if adj.ndim==1:
            if self.vID == 1:
                return gtool.normal_vec_2d(self._coords(adj[0]), self._coords(adj[1]))
            elif self.vID == 2:
                return gtool.normal_vec(self._coords(adj[0]),self._coords(adj[1]),self._coords(adj[2]))
            return
        v0 = np.array([adj[i][0] for i in range (adj.shape[0])])
        v1 = np.array([adj[i][1] for i in range (adj.shape[0])])
        if self.vID == 1:
            return gtool.normal_vec_2d(self._coords(v0),self._coords(v1))

            #edges_adj = self.connectivities[range_vec]
            #centers  = 0.5* (self._coords(edges_adj[:,0]) + self._coords(edges_adj[:,1]))
        elif self.vID == 2:
            v2 = np.array([adj[i][2] for i in range (adj.shape[0])])
            return  gtool.normal_vec(self._coords(v0),self._coords(v1),self._coords(v2))

    def _connectivities(self,index):
        if isinstance(index, np.int64) or isinstance(index, int):
            el_handle = np.array([self.elements_handle[index]])
        elif not isinstance(index, np.ndarray) and index is not None:
            el_handle = self.elements_handle[index].get_array()
        else:
            el_handle = self.elements_handle.get_array(index)
        return self.mb.get_ord_connectivity(el_handle, tag_handle = self.tag_handle)

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

    def init_partition(self):
        config = self.read_config('msCoarse.yml')
        particionador_type = config["Partitioner Scheme"]
        specific_attributes = config["Coarsening"]
        if particionador_type != '0':
            if self.dim == 3:
                partition = MoabVariable(self.core,data_size=1,var_type= "volumes",  data_format="int", name_tag="Partition",
                                             data_density="sparse")
                name_function = "scheme" + particionador_type
                used_attributes = []
                used_attributes.append(specific_attributes[0]["nx"])
                used_attributes.append(specific_attributes[1]["ny"])
                used_attributes.append(specific_attributes[2]["nz"])
                [partition[:],coarse_center]  = getattr(algoritmo, name_function)(self.volumes.center[:],
                           len(self), self.rx, self.ry, self.rz,*used_attributes)
            elif self.dim == 2:
                partition = MoabVariable(self.core,data_size=1,var_type= "faces",  data_format="int", name_tag="Partition",
                                             data_density="sparse")
                name_function = "scheme" + particionador_type
                specific_attributes = config["Coarsening"]
                used_attributes = []
                used_attributes.append(specific_attributes[0]["nx"])
                used_attributes.append(specific_attributes[1]["ny"])
                [partition[:],coarse_center]  = getattr(algoritmo, name_function)(self.faces.center[:],
                           len(self), self.rx, self.ry, self.rz,*used_attributes)
            return partition

    def _classify_element(self, index):
        range_vec = self.create_range_vec(index)
        range = self.range_index(range_vec)
        type_list = np.array([self.mb.type_from_handle(el) for el in range])
        return type_list

    def range_index(self, vec_index, flag_nodes=False):
        if not flag_nodes:
            range_handle = self.elements_handle
        else:
            range_handle = self.nodes
        if vec_index.dtype == "bool":
            vec = np.where(vec_index)[0]
        else:
            vec = vec_index.astype("uint")
        handles = np.asarray(range_handle)[vec.astype("uint")].astype("uint")
        return handles
        # return rng.Range(handles)

    def __str__(self):
        string = "{0} object \n Total of {1} {0} \n {2}  boundary {0} \n {3} internal {0}".format(self.entity_type,
            len(self.elements_handle), len(self.boundary_elements), len(self.internal_elements))
        return string

    def __len__(self):
        return len(self.elements_handle)

    def __call__(self):
        return self.all

    def read(self, handle):
        return self.mb.tag_get_data(self.tag_handle, handle, flat = True).astype(np.int64)

    @property
    def all_flagged_elements(self):
        return np.array(  list(self.flag.values())).astype(int)

    @property
    def all_flags(self):
        return np.array(list(self.flag.keys())).astype(int)

    @property
    def all(self):
        return self.read(self.elements_handle)

    @property
    def boundary(self):
        return self.read(self.boundary_elements.get_array())

    @property
    def internal(self):
        return self.read(self.internal_elements.get_array())


class MoabVariable(object):
    def __init__(self, core, name_tag, var_type="volumes", data_size=1, data_format="float", data_density="sparse",
                 entity_index=None):
        # pdb.set_trace()
        dic_dtf = {'float': types.MB_TYPE_DOUBLE, 'int': types.MB_TYPE_INTEGER, 'bool': types.MB_TYPE_BIT}
        dic_dst = {'dense': types.MB_TAG_DENSE, 'sparse': types.MB_TAG_SPARSE, 'bit': types.MB_TAG_BIT}
        dic_elem = {'nodes': core.all_nodes, 'edges': core.all_edges, 'faces': core.all_faces, 'volumes': core.all_volumes}
        self.mb = core.mb
        self.var_type = var_type
        self.data_format = data_format
        self.data_size = data_size
        self.data_density = data_density
        self.name_tag = name_tag
        self.custom = False
        if var_type in dic_elem:
            self.elements_handle = dic_elem[var_type]
        if entity_index is not None:
            self.elements_handle = self.range_index(entity_index)
            self.custom = True
        if data_density not in dic_dst:
            print("Please define a valid tag type")
        if data_format not in dic_dtf:
            print("Please define a valid data format")
        self.tag_handle = self.mb.tag_get_handle(name_tag, data_size, dic_dtf[data_format], dic_dst[data_density], True)
        print("Component class {0} successfully intialized".format(self.name_tag))

    def __call__(self):
        return self.mb.tag_get_data(self.tag_handle, self.elements_handle)

    def __setitem__(self, index, data):
        if isinstance(index, np.int64) or isinstance(index, int):
            el_handle = np.array([self.elements_handle[index]])
        elif not isinstance(index, np.ndarray) and index is not None:
            el_handle = self.elements_handle[index].get_array()
        else:
            el_handle = self.elements_handle.get_array(index)
        if isinstance(data, int) or isinstance(data, float) or isinstance(data, bool):
            data = data * np.ones((el_handle.size, self.data_size)).astype(self.data_format)
        elif (isinstance(data, np.ndarray)) and (len(data) == self.data_size):
            data = np.tile(data, (el_handle.size, 1)).astype(self.data_format)
        elif isinstance(data, list) & (len(data) == self.data_size):
            data = np.array(data)
            data = np.tile(data, (el_handle.size, 1)).astype(self.data_format)
        self.mb.tag_set_data(self.tag_handle, el_handle, data)

    def __getitem__(self, index):
        if isinstance(index, np.int64) or isinstance(index, int):
            el_handle = np.array([self.elements_handle[index]])
        elif not isinstance(index, np.ndarray) and index is not None:
            el_handle = self.elements_handle[index].get_array()
        else:
            el_handle = self.elements_handle.get_array(index)
        return self.mb.tag_get_data(self.tag_handle, el_handle)

    def __str__(self):
        string = "{0} variable: {1} based - Type: {2} - Length: {3} - Data Type: {4}"\
            .format(self.name_tag.capitalize(), self.var_type.capitalize(), self.data_format.capitalize(),
                    self.data_size, self.data_density.capitalize())
        if self.custom:
            string = string + " - Custom variable"
        return string

    def __len__(self):
        return len(self.elements_handle)
