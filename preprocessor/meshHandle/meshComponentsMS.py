"""
Generator of multiscale mesh entities and tags
"""
import numpy as np
from . meshComponents import MoabVariable, MeshEntities
from pymoab import types, rng
#mport pdb


class GetItem(object):
    def __init__(self, adj):
        self.fun = adj

    def __call__(self, item):
        return self.fun(item)

    def __getitem__(self, item):
        return self.fun(item)


class MeshEntitiesMS(MeshEntities):
    def __init__(self, core, entity_type):
        # pdb.set_trace()
        super().__init__(core, entity_type)
        # print(core)
        # print(entity_type)
    def enhance(self,i, general):
        self.coarse_neighbors_dic = {}
        if self.vID == 0:
            self.coarse_neighbors_dic = { key[1] :value for key, value in general.nodes_neighbors.items() if key[0] == i}
        elif self.vID == 1:
            self.coarse_neighbors_dic = {key[1]:value for key, value in general.edges_neighbors.items() if key[0] == i}
        elif self.vID == 2:
            self.coarse_neighbors_dic = {key[1]:value for key, value in general.faces_neighbors.items() if key[0] == i}
        elif self.vID == 3:
            self.coarse_neighbors_dic = {key[1]:value for key, value in general.volumes_neighbors.items() if key[0] == i}
        self.coarse_neighbors =  np.array([key for key, value in self.coarse_neighbors_dic.items() if not value.empty()])
        self.all_coarse_neighbors_range= rng.Range()
        for el in self.coarse_neighbors_dic.values():
            self.all_coarse_neighbors_range = rng.unite(self.all_coarse_neighbors_range,el)
        self.elements_in_coarse_neighborhood = GetItem(self._elements_in_coarse_neighborhood)

    def _elements_in_coarse_neighborhood(self,x):
        handle = self.coarse_neighbors_dic[x]
        return self.read(handle)

    @property
    def all_elements_in_coarse_neighborhood(self):
            return self.read(self.all_coarse_neighbors_range)


class MoabVariableMS(MoabVariable):
    def __init__(self, core, name_tag, var_type="volumes", data_size=1, data_format="float", data_density="sparse",
                 entity_index=None, level = 0, coarse_num = 0):
        self.mb = core.mb
        self.var_type = var_type
        self.data_format = data_format
        self.data_size = data_size
        self.data_density = data_density
        self.name_tag = name_tag
        self.custom = False
        if var_type == "nodes":
            self.elements_handle = core.all_nodes
        elif var_type == "edges":
            self.elements_handle = core.all_edges
        elif var_type == "faces":
            self.elements_handle = core.all_faces
        elif var_type == "volumes":
            self.elements_handle = core.all_volumes
        if entity_index is not None:
            self.elements_handle = self.range_index(entity_index)
            self.custom = True
        if data_density == "dense":
            data_density = types.MB_TAG_DENSE
        elif data_density == "sparse":
            data_density = types.MB_TAG_SPARSE
        elif data_density == "bit":
            data_density = types.MB_TAG_BIT
        else:
            print("Please define a valid tag type")
        if data_format == 'float':
            data_format = types.MB_TYPE_DOUBLE
        elif data_format == "int":
            data_format = types.MB_TYPE_INTEGER
        elif data_format == "bool":
            data_format = types.MB_TYPE_BIT
        self.level = level
        self.coarse_num = coarse_num
        self.name_tag = self.name_tag  + "-L" + str(self.level) + "-" + str(self.coarse_num)

        self.tag_handle = self.mb.tag_get_handle(self.name_tag, data_size, data_format, data_density, True)

        print("Component class {0} successfully intialized".format(self.name_tag))


class MultiscaleMeshEntities(object):
    def __init__(self,father_core,coarse_list):
        self.mb = father_core.mb
        self.num_coarse = len(coarse_list)
        self.find_coarse_neighbours(coarse_list)
        self.num = {"nodes": 0, "node": 0, "edges": 1, "edge": 1, "faces": 2, "face": 2, "volumes": 3, "volume": 3,
                             0: 0, 1: 1, 2: 2, 3: 3}
        # self.local_tag = coarse_list[0].core.handleDic[coarse_list[0].core.id_name]
        self.local_tag =  [el.core.handleDic["LOCAL_ID_L" + str(el.core.level) + "-" + str(el.core.coarse_num)] for el in coarse_list]
        self.global_tag = father_core.handleDic["GLOBAL_ID"]
        self.all_volumes = father_core.all_volumes
        self.all_faces = father_core.all_faces
        self.all_edges = father_core.all_edges
        self.all_nodes = father_core.all_nodes
        self.create_coarse_connectivities()

        #self.num = {"nodes": self.nodes_neighbors, "node": self.nodes_neighbors, "edges": self.edges_neighbors  , "edge": self.edges_neighbors , "faces": self.faces_neighbors, "face": self.faces_neighbors ,
        #                        "volumes": self.volumes_neighbors, "volume": self.volumes_neighbors, 0: self.nodes_neighbors , 1: self.edges_neighbors, 2: self.faces_neighbors , 3: self.all_volumes}

        # self.id_name = "LOCAL_ID_L" + str(self.level) + "-" + str(self.coarse_num)

        #pass

    def create_coarse_connectivities(self):
        self.connectivities = np.zeros((self.num_coarse,self.num_coarse))
        for x in range(self.num_coarse):
            for y in range(x+1,self.num_coarse):
                if not self.nodes_neighbors[x,y].empty():
                    self.connectivities[x,y]  =  True
                    self.connectivities[y,x]  =  True

    def find_coarse_neighbours(self,coarse_list):
        self.nodes_neighbors  = {}
        self.edges_neighbors  = {}
        self.faces_neighbors  = {}
        self.volumes_neighbors  = {}
        self.all_nodes_neighbors = rng.Range()
        self.all_edges_neighbors = rng.Range()
        self.all_faces_neighbors = rng.Range()
        self.all_volumes_neighbors = rng.Range()
        for x in range(self.num_coarse):
            for y in range(x+1,self.num_coarse):
                self.nodes_neighbors[x,y] = rng.intersect(coarse_list[x].core.boundary_nodes, coarse_list[y].core.boundary_nodes)
                self.nodes_neighbors[y,x] = self.nodes_neighbors[x,y]
                temp = self.nodes_neighbors[x,y]
                [self.all_nodes_neighbors.insert(e) for e in temp]

                self.edges_neighbors[x,y] = rng.intersect(coarse_list[x].core.boundary_edges, coarse_list[y].core.boundary_edges)
                self.edges_neighbors[y,x] = self.edges_neighbors[x,y]
                temp = self.edges_neighbors[x,y]
                [self.all_edges_neighbors.insert(e) for e in temp]

                self.faces_neighbors[x,y] = rng.intersect(coarse_list[x].core.boundary_faces, coarse_list[y].core.boundary_faces)
                self.faces_neighbors[y,x] = self.faces_neighbors[x,y]
                temp = self.faces_neighbors[x,y]
                [self.all_faces_neighbors.insert(e) for e in temp]

                self.volumes_neighbors[x,y] = rng.intersect(coarse_list[x].core.boundary_volumes, coarse_list[y].core.boundary_volumes)
                self.volumes_neighbors[y,x] = self.volumes_neighbors[x,y]
                temp = self.volumes_neighbors[x,y]
                [self.all_volumes_neighbors.insert(e) for e in temp]

    def global_to_local_id(self,vec_range,element, target ):
        flag = self.num[element]
        vec = self.create_range_vec(vec_range)
        if flag == 0:
            handle = self.range_index(vec, self.all_nodes)
        elif flag == 1:
            handle = self.range_index(vec, self.all_edges)
        elif flag == 2:
            handle = self.range_index(vec, self.all_faces)
        elif flag == 3:
            handle = self.range_index(vec, self.all_volumes)
        return self.mb.tag_get_data(self.local_tag[target],handle)

    def coarse_neighbours(self, x,y, element):
          # return self.read_data(self.global_tag, range_el = self.num[element])
          flag = self.num[element]
          if flag == 0:
              return self.mb.tag_get_data(self.global_tag, self.nodes_neighbors[x,y])
          elif flag == 1:
              return self.mb.tag_get_data(self.global_tag, self.edges_neighbors[x,y])
          elif flag == 2:
              return self.mb.tag_get_data(self.global_tag, self.faces_neighbors[x,y])
          elif flag == 3:
              return self.mb.tag_get_data(self.global_tag, self.volumes_neighbors[x,y])

    @property
    def all_neighbors_nodes(self):
        return self.mb.tag_get_data(self.global_tag, self.all_nodes_neighbors)

    @property
    def all_neighbors_edges(self):
        return self.mb.tag_get_data(self.global_tag, self.all_edges_neighbors)

    @property
    def all_neighbors_faces(self):
        return self.mb.tag_get_data(self.global_tag, self.all_faces_neighbors)

    @property
    def all_neighbors_volumes(self):
        return self.mb.tag_get_data(self.global_tag, self.all_volumes_neighbors)

    def create_range_vec(self, index):
        range_vec = None
        if isinstance(index, int):
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
