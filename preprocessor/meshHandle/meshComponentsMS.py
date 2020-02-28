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
        super().__init__(core, entity_type)
        self.father_handle = core.handleDic[core.father_id_name]
        self.global_id = GetItem(self._global_id)
        self.father_id = GetItem(self._father_id)

    def _global_id(self, index):
        el_handle = self.get_range_array(index)
        return self.mb.tag_get_data(self.global_handle, el_handle, flat = True).astype(np.int64)

    def _father_id(self, index):
        el_handle = self.get_range_array(index)
        return self.mb.tag_get_data(self.father_handle, el_handle, flat = True).astype(np.int64)

    def enhance(self,i, general):
        self._coarse_neighbors_dic = {}
        if self.vID == 0:
            index = general._nodes_neighbors[i, np.where(general.connectivities[i, :, 0])[0]].astype("uint64")
            self._coarse_neighbors_dic = {x: general._nodes[y] for x,y in zip(np.where(general.connectivities[i, :, 0])[0], index)}
        elif self.vID == 1:
            index = general._edges_neighbors[i, np.where(general.connectivities[i, :, 1])[0]].astype("uint64")
            self._coarse_neighbors_dic = {x: general._edges[y] for x,y in zip(np.where(general.connectivities[i, :, 1])[0], index)}
        elif self.vID == 2:
            index = general._faces_neighbors[i, np.where(general.connectivities[i, :, 2])[0]].astype("uint64")
            self._coarse_neighbors_dic = {x: general._faces[y] for x,y in zip(np.where(general.connectivities[i, :, 2])[0], index)}
        if self.vID < 3:
            self.coarse_neighbors = np.where(general.connectivities[i, :, self.vID])[0].astype("uint64")
            self.is_on_father_boundary = general.connectivities[i, -1, self.vID]
        self.neighborhood = GetItem(self._elements_in_coarse_neighborhood)

    def _elements_in_coarse_neighborhood(self,x):
        handle = self._coarse_neighbors_dic[x]
        return self.read(handle)

    @property
    def father_boundary(self):
        if self.is_on_father_boundary:
            handle = self.coarse_neighbors_dic[max(self.coarse_neighbors_dic.keys())]
            return self.read(handle)

    @property
    def all_coarse_neighbors(self):
        trange = rng.Range()
        for el in self.coarse_neighbors_dic.values():
            trange = rng.unite(trange, el)
        if self.is_on_father_boundary:
            trange = rng.subtract(trange, self.coarse_neighbors_dic[max(self.coarse_neighbors_dic.keys())])
        return self.read(trange)


class MoabVariableMS(MoabVariable):
    def __init__(self, core, name_tag, var_type="volumes", data_size=1, data_format="float", data_density="sparse",
                 entity_index=None, level=0, coarse_num=0, create = True):

        self.mb = core.mb
        self.var_type = var_type
        self.data_format = data_format
        self.data_size = data_size
        self.data_density = data_density
        self.name_tag = name_tag
        if var_type == "nodes":
            self.elements_handle = core.all_nodes
        elif var_type == "edges":
            self.elements_handle = core.all_edges
        elif var_type == "faces":
            self.elements_handle = core.all_faces
        elif var_type == "volumes":
            self.elements_handle = core.all_volumes
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
        #import pdb; pdb.set_trace()
        #"-L" + str(self.level) + "-" + str(self.coarse_num)

        self.storage = 'moab'
        self.moab_updated = True

        if create:
            name = core.id_name
            name = name[(name.find("ID") + 3):]
            self.name_tag = self.name_tag  + name
            self.tag_handle = self.mb.tag_get_handle(self.name_tag, data_size, data_format, data_density, True, 0)
        else:
            self.tag_handle = self.mb.tag_get_handle(self.name_tag)
        print("Component class {0} successfully intialized".format(self.name_tag))
