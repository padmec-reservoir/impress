"""
Generator of mesh entities and tags
"""

import numpy as np
from pymoab import types, rng, topo_util
from ..geoUtil import geoTools as gtool
from .get_item import GetItem

class MeshEntities(object):
    def __init__(self, core, entity_type):
        # Main MOAB structures.
        self.mb = core.mb
        self.mtu = core.mtu
        self.meshset = core.root_set

        self.entity_type = entity_type

        self.nodes = core.all_nodes

        # Map from string description to int.
        self.entity_str_to_num = {"nodes": 0, "node": 0, 
                    "edges": 1, "edge": 1, 
                    "faces": 2, "face": 2, 
                    "volumes": 3, "volume": 3,
                    0: 0, 1: 1, 2: 2, 3: 3}

        # List of elements of all dimensions.
        self.list_all = [core.all_nodes, core.all_edges, core.all_faces, core.all_volumes]

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
            self.id_name = (self.father_id_name + str("L") + str(self.level) + "-"
                + str(core.coarse_num))

        # Set mesh entities.
        self.vID = self.entity_str_to_num[entity_type]
        if self.vID == 0:
            self.elements_handle = core.all_nodes
            self.internal_range = core.internal_nodes
            self.boundary_range = core.boundary_nodes
        elif self.vID == 1:
            self.elements_handle = core.all_edges
            self.internal_range = core.internal_edges
            self.boundary_range = core.boundary_edges
        elif self.vID == 2:
            self.elements_handle = core.all_faces
            self.internal_range = core.internal_faces
            self.boundary_range = core.boundary_faces
        elif self.vID == 3:
            self.elements_handle = core.all_volumes
            self.internal_range = core.internal_volumes
            self.boundary_range = core.boundary_volumes
        else:
            raise ValueError("Invalid dimension for mesh entities.")

        self.tag_handle = core.handleDic[self.id_name]
        self.global_handle = core.handleDic['GLOBAL_ID']
        self.father_handle = core.handleDic[self.father_id_name]

        # Set mesh entities properties.
        self.all_elements = GetItem(self.get_all)
        self.internal_elements = GetItem(self.get_internal)
        self.boundary_elements = GetItem(self.get_boundary)
        if self.vID == 0:
            self.adjacencies = GetItem(self._adjacencies_for_nodes)
            self.coords =  GetItem(self._coords)
        else:
            self.adjacencies = GetItem(self._adjacencies)
            self.connectivities = GetItem(self._connectivities)
        self.classify_element = GetItem(self._classify_element)
        self.center = GetItem(self._center)
        self.area = GetItem(self._area)
        self.volume = GetItem(self._volume)
        if (self.vID == 1) & (core.dimension == 2):
            self.normal = GetItem(self._normal)
        elif (self.vID == 2) & (core.dimension == 3):
            self.normal = GetItem(self._normal)
        
        # Initialize specific flag dic according to the type of the object created.
        self.flag = {key: self.read(value[self.vID]) for key, value in core.flag_dic.items()
                     if value[self.vID].empty() is not True}

    def __str__(self):
        string = "{0} object \n Total of {1} {0} \n {2}  boundary {0} \n {3} internal {0}".format(self.entity_type,
            len(self.elements_handle), len(self.boundary_elements), len(self.internal_elements))
        return string

    def __len__(self):
        return len(self.elements_handle)

    def __call__(self):
        return self.all

    def format_entities(self, input_tuple, input_size, tag = None):
        entities, sizes, jagged = input_tuple
        if tag is not None:
            entities = self.mb.tag_get_data(tag, entities, flat = True).astype(np.int64)
        if jagged:
            entities = np.delete(np.array(np.split(entities, sizes), dtype=object), -1)
        else:
            entities = entities.reshape(-1, sizes)
        
        if input_size == 1:
            entities = entities[0]
        
        return entities
    
    def bridge_adjacencies(self, index, interface, target, ordering_inst=None):
        el_handle = self.get_range_array(index)
        intersect_ent = None

        if self.level > 0:
            intersect_ent = self.list_all[self.entity_str_to_num[target]]
        
        result_tuple = self.mtu.get_ord_bridge_adjacencies(el_handle, 
            self.entity_str_to_num[interface], self.entity_str_to_num[target], 
            intersect_ent, self.level)

        entities_array = self.format_entities(result_tuple, el_handle.size, self.tag_handle)

        if ordering_inst != None:
            entities_array = ordering_inst.sort_elements(entities_array, self.center[index])

        return entities_array

    def _coords(self, index):
        el_handle = self.get_range_array(index, self.nodes)
        if len(el_handle) == 1:
            return self.mb.get_coords(el_handle)
        return np.reshape(self.mb.get_coords(el_handle),(-1,3))

    def _adjacencies_for_nodes(self, index):
        return index

    def _adjacencies(self, index, flag_nodes=False, dim_tag=None):
        if dim_tag is None:
            if not flag_nodes:
                dim_tag = self.vID - 1
            else:
                dim_tag = 0
        el_handle = self.get_range_array(index)
        result_tuple = self.mb.get_ord_adjacencies(el_handle, dim_tag)
        return self.format_entities(result_tuple, el_handle.size, self.tag_handle)

    def _center(self, index):
        el_handle = self.get_range_array(index)
        return self.mtu.get_ord_average_position(el_handle)

    def _normal(self, index):
        adj = self.connectivities[index]
        normal = None

        if adj.ndim == 1 and self.vID == 1:
            normal = gtool.normal_vec_2d(self._coords(adj[0]).reshape(1,3), 
                                        self._coords(adj[1]).reshape(1,3))[0]
        elif adj.ndim == 1 and self.vID == 2:
            normal = gtool.normal_vec(self._coords(adj[0]), 
                                    self._coords(adj[1]), 
                                    self._coords(adj[2]))
        else:
            v0 = np.array([adj[i][0] for i in range (adj.shape[0])])
            v1 = np.array([adj[i][1] for i in range (adj.shape[0])])
            if self.vID == 1:
                normal = gtool.normal_vec_2d(self._coords(v0),self._coords(v1))
            elif self.vID == 2:
                v2 = np.array([adj[i][2] for i in range (adj.shape[0])])
                normal = gtool.normal_vec(self._coords(v0), self._coords(v1),
                                        self._coords(v2))
        
        return normal

    def _area(self, index):
        if self.entity_type != "faces":
            raise ValueError("Cannot compute area. Entity is not a face.")

        el_handle = self.get_range_array(index)[0]
        area = gtool.polygon_area(self.mb, el_handle)

        return area
    
    def _volume(self, index):
        if self.entity_type != "volumes":
            raise ValueError("Cannot compute volume. Entity is not a volume.")

        el_handle = self.get_range_array(index)[0]
        center = self._center(index)
        volume = gtool.polyhedron_volume(self.mb, el_handle, center)

        return volume

    def _connectivities(self,index):
        el_handle = self.get_range_array(index)
        result_tuple = self.mb.get_ord_connectivity(el_handle)
        return self.format_entities(result_tuple, el_handle.size, self.tag_handle)

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
            start, stop, step = index.start, index.stop, index.step
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

    def _classify_element(self, index):
        range_vec = self.create_range_vec(index)
        a_range = self.range_index(range_vec)
        type_list = np.array([self.mb.type_from_handle(el) for el in a_range])
        return type_list

    def load_array(self, array=None):
        if array == None:
            self.all_elements = self.all_elements[:]
            self.internal_elements = self.internal_elements[:]
            self.boundary_elements = self.boundary_elements[:]
        elif array == 'all':
            self.all_elements = self.all_elements[:]
        elif array == 'internal':
            self.internal_elements = self.internal_elements[:]
        elif array == 'boundary':
            self.boundary_elements = self.boundary_elements[:]

    def unload_array(self, array=None):
        if array == None:
            self.all_elements = GetItem(self.get_all)
            self.internal_elements = GetItem(self.get_internal)
            self.boundary_elements = GetItem(self.get_boundary)
        elif array == 'all':
            self.all_elements = GetItem(self.get_all)
        elif array == 'internal':
            self.internal_elements = GetItem(self.get_internal)
        elif array == 'boundary':
            self.boundary_elements = GetItem(self.get_boundary)

    def get_range_array(self, index, search_range=None):
        el_handle = None
        if search_range == None:
            search_range = self.elements_handle
        if isinstance(index, np.int64) or isinstance(index, int):
            return np.array([search_range[index]])
        if not isinstance(index, np.ndarray) and index is not None:
            el_handle = search_range[index].get_array()
        else:
            el_handle = search_range.get_array(index)
        return el_handle

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

    def read(self, handle):
        return self.mb.tag_get_data(self.tag_handle, handle, flat = True).astype(np.int64)

    @property
    def all_flagged_elements(self):
        return np.array(list(self.flag.values())).astype(int)

    @property
    def all_flags(self):
        return np.array(list(self.flag.keys())).astype(int)

    @property
    def all(self):
        return self.read(self.elements_handle.get_array())

    @property
    def boundary(self):
        return self.read(self.boundary_range.get_array())

    @property
    def internal(self):
        return self.read(self.internal_range.get_array())

    def get_all(self, index):
        if self.level == 0 and isinstance(index, np.ndarray):
            ret = index.astype(np.int64)
        elif self.level == 0 and isinstance(index, slice):
            if index.start == None and index.stop == None:
                ret = np.arange(self.elements_handle.size())
            elif index.start == None:
                ret = np.arange(0, min(index.stop, self.elements_handle.size()), 
                                index.step, dtype = np.int64)
            elif index.stop == None:
                ret = np.arange(index.start, self.elements_handle.size(), 
                                index.step, dtype = np.int64)
            else:
                ret = np.arange(index.start, min(index.stop, self.elements_handle.size()), 
                                index.step, dtype = np.int64)
        else:
            ret = self.read(self.elements_handle[index])
        
        if len(ret) == 1:
            ret = ret[0]
        
        return ret

    def get_boundary(self, index):
        el_range = self.get_range_array(index, self.boundary_range)
        ret = self.read(el_range)

        if len(ret) == 1:
            ret = ret[0]
        
        return ret

    def get_internal(self, index):
        el_range = self.get_range_array(index, self.internal_range)
        ret = self.read(el_range)

        if len(ret) == 1:
            ret = ret[0]
        
        return ret

class MoabVariable(object):
    def __init__(self, core, name_tag, var_type="volumes", 
                data_size=1, data_format="float", 
                data_density="sparse", entity_index=None, 
                create = True):
        dic_data_format = {'float': types.MB_TYPE_DOUBLE, 
                            'int': types.MB_TYPE_INTEGER, 
                            'bool': types.MB_TYPE_BIT}
        dic_density = {'dense': types.MB_TAG_DENSE, 
                        'sparse': types.MB_TAG_SPARSE, 
                        'bit': types.MB_TAG_BIT}
        dic_elem = {'node': core.all_nodes, 
                    'nodes': core.all_nodes, 
                    'edge': core.all_edges, 
                    'edges': core.all_edges, 
                    'face': core.all_faces, 
                    'faces': core.all_faces, 
                    'volume': core.all_volumes, 
                    'volumes': core.all_volumes}
        
        self.data = None
        self.mb = core.mb
        self.var_type = var_type
        self.data_format = data_format
        self.data_size = data_size
        self.data_density = data_density
        self.name_tag = name_tag

        if var_type in dic_elem:
            self.elements_handle = dic_elem[var_type]
        else:
            raise ValueError('var_type must be "nodes", "edges", "faces" or "volumes"')

        if data_density not in dic_density:
            raise ValueError('Data density must be "dense", "sparse" or "bit".')
        
        if data_format not in dic_data_format:
            raise ValueError('Data format must be "float", "int" or "bool".')
        
        if create == True:
            self.tag_handle = self.mb.tag_get_handle(name_tag, data_size, 
                                                    dic_data_format[data_format], 
                                                    dic_density[data_density], 
                                                    True, 0)
        else:
            self.tag_handle = self.mb.tag_get_handle(name_tag)
        self.storage = 'moab'
        self.moab_updated = True
        print("Component class {0} successfully intialized".format(self.name_tag))
    
    def __call__(self):
        ret_value = None
        if self.storage == 'moab':
            ret_value = self.mb.tag_get_data(self.tag_handle, self.elements_handle)
        else:
            ret_value = self.data
        return ret_value

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

        if self.storage == 'moab':
            self.mb.tag_set_data(self.tag_handle, el_handle, data)

        if self.storage == 'numpy':
            if self.data_size == 1:
                data = data.flatten()
            self.data[index] = data
            self.moab_updated = False

    def __getitem__(self, index=None):
        item = None

        if self.storage == 'moab':
            if isinstance(index, np.int64) or isinstance(index, int):
                el_handle = np.array([self.elements_handle[index]])
            elif not isinstance(index, np.ndarray) and index is not None:
                el_handle = self.elements_handle[index].get_array()
            else:
                el_handle = self.elements_handle.get_array(index)
            item = self.mb.tag_get_data(self.tag_handle, el_handle)
        else:
            if index is not None:
                item = self.data[index]
            else:
                item = self.data
        
        return item
    
    def __str__(self):
        string = "{0} variable: {1} based - Type: {2} - Length: {3} - Data Type: {4}"\
            .format(self.name_tag.capitalize(), self.var_type.capitalize(), 
                    self.data_format.capitalize(), self.data_size, 
                    self.data_density.capitalize())
        return string

    def __len__(self):
        return len(self.elements_handle)

    def to_numpy(self):
        if self.storage == 'numpy':
            raise UserWarning("Variable is already on numpy")
        self.storage = 'numpy'
        self.data = self.mb.tag_get_data(self.tag_handle, self.elements_handle)
        if self.data_size == 1:
            self.data = self.data.flatten()

    def to_moab(self):
        if self.storage == 'moab':
            raise UserWarning("Variable is already on moab")
        self.storage = 'moab'
        self.mb.tag_set_data(self.tag_handle, self.elements_handle, self.data)
        self.data = None
        self.moab_updated = True

    def update_moab(self):
        self.mb.tag_set_data(self.tag_handle, self.elements_handle, self.data)
        self.moab_updated = True
