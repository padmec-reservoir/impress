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
        name = core.id_name
        name = name[(name.find("ID") + 3):]
        self.name_tag = self.name_tag  + name
        #import pdb; pdb.set_trace()
        print(self.name_tag)
        print(data_size)
        print(data_format)
        print(data_density)
        #"-L" + str(self.level) + "-" + str(self.coarse_num)
        self.tag_handle = self.mb.tag_get_handle(self.name_tag, data_size, data_format, data_density, True)
        print("Component class {0} successfully intialized".format(self.name_tag))
