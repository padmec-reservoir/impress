"""
Use of Pymoab methods to manage the multiscale mesh
"""
#import pdb
from . corePymoab import CoreMoab
from pymoab import core, types, rng
import numpy as np

class MsCoreMoab(CoreMoab):
    def __init__(self, father_core, num, coarse_vec):
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
            self.all_faces = self.access_handle(self.all_volumes)
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
        self.id_name = "LOCAL_ID_L" + str(self.level) + "-" + str(self.coarse_num)
        self.init_id()
        # all_entities = self.mb.get_entities_by_handle(self.root_set)
        # print(all_entities)
        self.flag_dic = {key:[rng.intersect(all_entities,el) for el in value] for (key, value) in father_core.flag_dic.items()}

    def skinner_operation(self):
        #skin = sk.Skinner(self.mb)
        # print("Entering skinner test")

        if self.dimension == 3:
            # faces_on_skin_handles = skin.find_skin(self.root_set, self.all_volumes)
            # faces_on_skin_handles = self.bridge_adjacencies(self.all_faces,self.dimension)

            # pdb.set_trace()
            faces_on_skin_handles = self.bridge_adjacencies(self.all_faces,self.dimension)
            # pdb.set_trace()

            # pdb.set_trace()
            edges_on_skin_handles = self.access_handle(faces_on_skin_handles)

            nodes_on_skin_handles = self.access_handle(edges_on_skin_handles)


            nodes_in_volumes = ([self.mb.get_adjacencies(el_handle,0) for el_handle in self.all_volumes])
            # check_volumes = np.asarray([(rng.intersect(el_handle,nodes_on_skin_handles))  for el_handle in nodes_in_volumes])
            check_volumes = ([(rng.intersect(el_handle,nodes_on_skin_handles))  for el_handle in nodes_in_volumes])


            external_volumes_index = np.array([el_handle.empty() for el_handle in check_volumes]).astype(bool)
            volumes_on_skin_handles = self.range_index(np.bitwise_not(external_volumes_index),self.all_volumes)


        elif self.dimension == 2:

            edges_on_skin_handles = self.bridge_adjacencies(self.all_edges,self.dimension)


            # edges_on_skin_handles = skin.find_skin(self.root_set, self.all_faces[:])
            nodes_on_skin_handles = self.access_handle(edges_on_skin_handles)
            nodes_in_faces = ([self.mb.get_adjacencies(el_handle,0) for el_handle in self.all_faces])
            check_faces= ([rng.intersect(el_handle,nodes_on_skin_handles) for el_handle in nodes_in_faces])
            external_faces_index = np.array([el_handle.empty() for el_handle in check_faces]).astype(bool)
            faces_on_skin_handles = self.range_index(np.bitwise_not(external_faces_index),self.all_faces)
            volumes_on_skin_handles = rng.Range()

        # print("Skinning Operation Successful")
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



#
#
# class MsCoreMoab:
#     def __init__(self, father_core, coarse_vec):
#         self.dimension = father_core.dimension
#         self.mb = father_core.mb
#         self.level = father_core.level + 1
#         self.father_root_set = father_core.root_set
#         self.root_set = self.mb.create_meshset(types.MESHSET_TRACK_OWNER)
#
#         self.mtu = father_core.mtu
#         #
#         self.handleDic = father_core.handleDic
#
#         self.all_volumes = self.range_index(coarse_vec, range_handle=father_core.all_volumes)
#         self.all_faces = self.access_handle(self.all_volumes)
#         self.all_edges = self.access_handle(self.all_faces)
#         self.all_nodes = rng.Range(self.mb.get_connectivity(self.all_volumes))
#
#         self.mb.add_entities(self.root_set, self.all_volumes)
#         self.mb.add_entities(self.root_set, self.all_faces)
#         self.mb.add_entities(self.root_set, self.all_edges)
#         self.mb.add_entities(self.root_set, self.all_nodes)
#         [self.boundary_nodes, self.boundary_edges, self.boundary_faces, self.boundary_volumes] = self.skinner_operation()
#         #
#         self.internal_nodes = rng.subtract(self.all_nodes, self.boundary_nodes)
#         self.internal_edges = rng.subtract(self.all_edges, self.boundary_edges)
#         self.internal_faces = rng.subtract(self.all_faces, self.boundary_faces)
#         self.internal_volumes = rng.subtract(self.all_volumes, self.boundary_volumes)
#         #
#         self.id_name = "LOCAL_ID_L" + str(self.level)
#         self.init_id()
#         all_entities = self.mb.get_entities_by_handle(self.root_set)
#         self.flag_dic = {key:[rng.intersect(all_entities,el) for el in value] for (key, value) in father_core.flag_dic.items()}
#
#
#
#
#
#
#
#         #self.mb.get_entities_by_type_and_tag(self.father_root_set, father_core.all_volumes,self.handleDic['GLOBAL_ID'],np.arange(len(father_core.all_volumes),dtype=int))
#
#         # self.root_set = self.mb.get_root_set()
#         # self.mtu = topo_util.MeshTopoUtil(self.mb)
#         # self.mb.load_file(mesh_file)
#         # self.all_volumes = self.mb.get_entities_by_dimension(0, 3)
#         # self.all_nodes = self.mb.get_entities_by_dimension(0, 0)
#         # self.mtu.construct_aentities(self.all_nodes)
#         # self.all_faces = self.mb.get_entities_by_dimension(0, 2)
#         # self.all_edges = self.mb.get_entities_by_dimension(0, 1)
#         # self.handleDic = {}
#         [self.boundary_nodes, self.boundary_edges, self.boundary_faces, self.boundary_volumes] = self.skinner_operation()
#         #
#         self.internal_nodes = rng.subtract(self.all_nodes, self.boundary_nodes)
#         self.internal_edges = rng.subtract(self.all_edges, self.boundary_edges)
#         self.internal_faces = rng.subtract(self.all_faces, self.boundary_faces)
#         self.internal_volumes = rng.subtract(self.all_volumes, self.boundary_volumes)
#         #
#         self.init_id()
#         all_entities = self.mb.get_entities_by_handle(self.root_set)
#          # rng.intersect
#         self.flag_dic = { key:[rng.intersect(all_entities,el) for el in value] for (key, value) in father_core.flag_dic.items()}
#         #
#         #
#         # for (key, value) in father_core.flag_dic.items():
#         #     a = [rng.intersect(all_entities,el) for el in value]
#         #     print(a)
#         #
#         #
#         # for (key, value) in father_core.flag_dic.items():
#         #     tmp_list = []
#         #     for el in value:
#         #         tmp_list.append(rng.intersect(all_entities,el))
#         #     print(tmp_list)
#             # tmp_list = []
#             # for el in value:
#             #     print(el)
#             #     tmp_list.append(rng.intersect(all_entities,el))
#             # self.flag_dic[key] = tmp_list
#
#         # elf.flag_dic = {key:[rng.intersect(all_entities,el)] for el in value for (key, value) in father_core.flag_dic.items()}
#
#         # for handle in values for key, value in father_core.flag_dic.items()]
#
#
#         #self.flag = {key: self.read(value[self.vID]) for key, value in core.flag_dic.items()
#         #             if value[self.vID].empty() is not True}
#
#         # self.check_integrity()
#         # self.create_id_visualization()
#         # # self.flag_dic = {}
#         # [self.flag_list, self.flag_dic] = self.read_flags()
#         #
#         # self.create_flag_visualization()
#         #
#         # # swtich on/off
#         #self.parallel_meshset = self.create_parallel_meshset()
#         #self.create_parallel_visualization()
#         # def init_entities(self, handle):
#         #     # read_data(self, name_tag, index_vec = np.array([]), range_el = None):
#         #
#         #     self.create_tag_handle('PARTITION', data_size, data_text = "int"):
#         #     m = self.read_data(self.handleDic['PARTITION'], range_el = handle)
#         #     pdb.set_trace()
#         #     pass
#
#     def init_id(self):
#         # delete previous IDs
#         # Gmesh standard counts from 1
#         # GLOBAL_ID_tag = self.mb.tag_get_handle(
#         #    "Global_ID", 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, False)
#
#         name_tag = "LOCAL_ID_L" + str(self.level)
#         # global_tag = self.mb.tag_get_handle(
#         #     name, 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, False)
#         # self.handleDic[name] = global_tag
#
#         self.create_tag_handle(name_tag, data_size = 1, data_text = "int", data_density = "sparse")
#         # create volume ids
#         self.set_data(name_tag, np.arange(len(self.all_volumes)))
#         # create face ids
#         self.set_data(name_tag, np.arange(len(self.all_faces)), range_el=self.all_faces)
#         # create edges ids
#         self.set_data(name_tag, np.arange(len(self.all_edges)), range_el=self.all_edges)
#         # create nodes ids
#         self.set_data(name_tag, np.arange(len(self.all_nodes)), range_el=self.all_nodes)
#
#     def skinner_operation(self):
#         skin = sk.Skinner(self.mb)
#         print("Entering skinner test")
#
#         if self.dimension == 3:
#             faces_on_skin_handles = skin.find_skin(self.root_set, self.all_volumes[:])
#             # pdb.set_trace()
#             edges_on_skin_handles = self.access_handle(faces_on_skin_handles)
#             nodes_on_skin_handles = self.access_handle(edges_on_skin_handles)
#
#             nodes_in_volumes = ([self.mb.get_adjacencies(el_handle,0) for el_handle in self.all_volumes])
#             check_volumes = np.asarray([rng.intersect(el_handle,nodes_on_skin_handles) for el_handle in nodes_in_volumes])
#             external_volumes_index = np.array([el_handle.empty() for el_handle in check_volumes]).astype(bool)
#             volumes_on_skin_handles = self.range_index(np.bitwise_not(external_volumes_index),self.all_volumes)
#
#
#         elif self.dimension == 2:
#             edges_on_skin_handles = skin.find_skin(self.root_set, self.all_faces[:])
#             nodes_on_skin_handles = self.access_handle(edges_on_skin_handles)
#             nodes_in_faces = ([self.mb.get_adjacencies(el_handle,0) for el_handle in self.all_faces])
#             check_faces= np.asarray([rng.intersect(el_handle,nodes_on_skin_handles) for el_handle in nodes_in_faces])
#             external_faces_index = np.array([el_handle.empty() for el_handle in check_faces]).astype(bool)
#             faces_on_skin_handles = self.range_index(np.bitwise_not(external_faces_index),self.all_faces)
#             volumes_on_skin_handles = rng.Range()
#
#         print("Skinning Operation Successful")
#         return [nodes_on_skin_handles, edges_on_skin_handles, faces_on_skin_handles, volumes_on_skin_handles]
#
#         #
#         # self.create_tag_handle("SKINPOINT",1, "int", data_density="sparse")
#         # self.create_tag_handle("SKINEDGES", 1, "int", data_density="sparse")
#         # self.create_tag_handle("SKINFACES", 1, "int", data_density="sparse")
#         # self.create_tag_handle("SKINVOLUMES", 1, "int", data_density="sparse")
#         #
#         # self.set_data("SKINPOINT", 200 * np.ones(len(nodes_on_skin_handles)).astype(int),range_el=nodes_on_skin_handles)
#         # self.set_data("SKINEDGES", 300 * np.ones(len(edges_on_skin_handles)).astype(int),range_el=edges_on_skin_handles)
#         # self.set_data("SKINFACES", 400 * np.ones(len(faces_on_skin_handles)).astype(int),range_el=faces_on_skin_handles)
#         # self.set_data("SKINVOLUMES", 800 * np.ones(len(volumes_on_skin_handles)).astype(int),range_el=volumes_on_skin_handles)
#         #
#         # return [nodes_on_skin_handles, edges_on_skin_handles, faces_on_skin_handles, volumes_on_skin_handles]
#
#     def check_integrity(self):
#         # check if the mesh contains
#         check_list = [len(self.all_nodes), len(self.all_edges), len(self.all_faces), len(self.all_volumes)]
#         list_words = ['Nodes', "Edges", "Faces", "Volumes"]
#         print("Checking mesh integrity:")
#         index = 0
#         if self.dimension == 2:
#             list_words = list_words[:-1]
#             check_list = check_list[:-1]
#         for entity in check_list:
#             if entity > 0:
#                 print(list_words[index] + " successfully imported")
#             else:
#                 print("------------------------------\nError creating \n" +
#                       list_words[index] + " was not imported")
#             index += 1
#
#     def create_parallel_meshset(self):
#         partition_volumes = []
#         parallel_tag = []
#         try:
#             parallel_tag = self.mb.tag_get_handle("PARALLEL_PARTITION")
#             flag = False
#         except:
#             print("Parallel Partition Tag not found \nAborting creating parallel partition entities.")
#             flag = True
#         if not flag:
#             print("Parallel Partition Tag detected \nCreating parallel partitions entities")
#             self.handleDic["PARALLEL_PARTITION"] = parallel_tag
#             parallel_sets = self.mb.get_entities_by_type_and_tag(
#                 0, types.MBENTITYSET, np.array(
#                 (parallel_tag,)), np.array((None,)))
#             self.create_tag_handle("PARALLEL", data_size=1, data_text="int")
#             # partition_volumes = []
#             for set_el in parallel_sets:
#                 num_tag = self.read_data("PARALLEL_PARTITION", range_el=set_el)[0, 0]
#                 list_entity = [self.mb.get_entities_by_dimension(set_el, 0), self.mb.get_entities_by_dimension(set_el, 1),
#                                self.mb.get_entities_by_dimension(set_el, 2), self.mb.get_entities_by_dimension(set_el, 3)]
#                 # print([num_tag, list_entity])
#                 partition_volumes.append(list_entity)
#         return partition_volumes
#
#     def create_parallel_visualization(self):
#         k = 0
#         for sets in self.parallel_meshset:
#             for dim in sets:
#                 if len(dim) != 0:
#                     self.set_data("PARALLEL", k * np.ones(len(dim)).astype(int), range_el=dim)
#             k += 1
#
#     def read_flags(self):
#         physical_tag = self.mb.tag_get_handle("MATERIAL_SET")
#         physical_sets = self.mb.get_entities_by_type_and_tag(
#             0, types.MBENTITYSET, np.array(
#             (physical_tag,)), np.array((None,)))
#         self.handleDic["MATERIAL_SET"] = physical_tag
#         flag_list = np.array([])
#         flag_dic = {}
#         for set in physical_sets:
#             bc_flag = self.read_data("MATERIAL_SET", range_el=set)[0, 0]
#             flag_list = np.append(flag_list, bc_flag)
#             list_entity = [self.mb.get_entities_by_dimension(set, 0), self.mb.get_entities_by_dimension(set, 1),
#                            self.mb.get_entities_by_dimension(set, 2), self.mb.get_entities_by_dimension(set, 3)]
#             flag_dic[bc_flag] = list_entity
#         return np.sort(flag_list), flag_dic
#
#     def create_flag_visualization(self):
#         self.create_tag_handle("FLAGS", 1, data_text="int")
#         self.create_tag_handle("MATERIAL", 1, data_text="int")
#         for k in self.flag_dic:
#             sets = self.flag_dic[k]
#             for dim in sets:
#                 if len(dim) != 0:
#                     if k > 100:
#                         self.set_data("FLAGS", k * np.ones(len(dim)).astype(int), range_el=dim)
#                     else:
#                         self.set_data("MATERIAL", k * np.ones(len(dim)).astype(int), range_el=dim)
#
#     def create_flag_visualization_alternative(self):
#         self.create_tag_handle("FLAGS-NODES", 1, data_text="int")
#         self.create_tag_handle("FLAGS-EDGES", 1, data_text="int")
#         self.create_tag_handle("FLAGS-FACES", 1, data_text="int")
#         self.create_tag_handle("FLAGS-VOUMES", 1, data_text="int")
#         for k in self.flag_dic:
#             sets = self.flag_dic[k]
#             for dim, ndim in zip(sets,range(4)):
#
#                 if len(dim) != 0:
#                     print([k, ndim])
#                     if k > 100:
#                         if ndim == 0:
#                             # print([dim, k])
#                             self.set_data("FLAGS-NODES", k * np.ones(len(dim)).astype(int), range_el=dim)
#                         elif ndim == 1:
#                             # print([dim, k])
#                             self.set_data("FLAGS-EDGES", k * np.ones(len(dim)).astype(int), range_el=dim)
#                         elif ndim ==2:
#                             # print([dim, k])
#                             self.set_data("FLAGS-FACES", k * np.ones(len(dim)).astype(int), range_el=dim)
#                         elif ndim == 3:
#                             # print([dim, k])
#                             self.set_data("FLAGS-VOLUMES", k * np.ones(len(dim)).astype(int), range_el=dim)
#                     # else:
#                     #     # self.set_data("MATERIAL", k * np.ones(len(dim)).astype(int), range_el=dim)
#                     #     pass
#
#     def create_id_visualization(self):
#         self.create_tag_handle("ID-NODES", 1, data_text="int")
#         self.create_tag_handle("ID-EDGES", 1, data_text="int")
#         self.create_tag_handle("ID-FACES", 1, data_text="int")
#         self.create_tag_handle("ID-VOLUMES", 1, data_text="int")
#
#         data_node = self.read_data("GLOBAL_ID", range_el=self.all_nodes)
#         data_edges = self.read_data("GLOBAL_ID", range_el=self.all_edges)
#         data_faces = self.read_data("GLOBAL_ID", range_el=self.all_faces)
#         data_volumes = self.read_data("GLOBAL_ID", range_el=self.all_volumes)
#
#         self.set_data("ID-NODES", data_node, range_el=self.all_nodes)
#         self.set_data("ID-EDGES", data_edges, range_el=self.all_edges)
#         self.set_data("ID-FACES", data_faces, range_el=self.all_faces)
#         self.set_data("ID-VOLUMES", data_volumes, range_el=self.all_volumes)
#
#     def access_meshset(self, handle):
#         # returns the entities contained inside a give meshset handle
#         # ie: for a meshset handle the entities inside are returned
#         temp_range = []
#         for el in range(self.dimension+1):
#             sub_el = (self.mb.get_entities_by_dimension(handle, el))
#             temp_range.append(sub_el)
#         temp_range.append(self.mb.get_entities_by_dimension(handle, 11))
#         return temp_range
#
#     def access_handle(self,handle):
#         # input: range of handles of different dimensions
#         # returns all entities with d-1 dimension the comprises the given range
#         # ie: for a volume, the faces, for a face the edges and for an edge the points.
#         #
#         vecdim = self.check_range_by_dimm(handle)
#         # pdb.set_trace()
#         all_adj = np.array([np.array(self.mb.get_adjacencies(el_handle, dim-1)) for dim, el_handle in zip(vecdim,handle)])
#         #unique_adj = np.unique(np.ma.concatenate(all_adj)).astype("uint64")
#         unique_adj = np.unique(np.concatenate(all_adj)).astype("uint64")
#         return rng.Range(unique_adj)
#
#     def create_tag_handle(self, name_tag, data_size, data_text = "float", data_density = "dense"):
#         if data_density == "dense":
#             data_density = types.MB_TAG_DENSE
#         elif data_density == "sparse":
#             data_density = types.MB_TAG_SPARSE
#         elif data_density == "bit":
#             data_density = types.MB_TAG_BIT
#         else:
#             print("Please define a valid tag type")
#         if data_text == 'float':
#             data_type = types.MB_TYPE_DOUBLE
#         elif data_text == "int":
#             data_type = types.MB_TYPE_INTEGER
#         elif data_text == "bool":
#             data_type = types.MB_TYPE_BIT
#         try:
#             handle = self.handleDic[name_tag]
#         except KeyError:
#             handle = self.mb.tag_get_handle(name_tag, data_size, data_type, data_density, True)
#             self.handleDic[name_tag] = handle
#
#     def read_data(self, name_tag, index_vec = np.array([]), range_el = None):
#         if range_el is None:
#             range_el = self.all_volumes
#         if index_vec.size > 0:
#             range_el = self.range_index(index_vec,range_el)
#         try:
#             handle_tag = self.handleDic[name_tag]
#             return self.mb.tag_get_data(handle_tag, range_el)
#         except KeyError:
#             print("Tag not found")
#
#     def range_index(self, vec_index, range_handle=None):
#         if range_handle is None:
#             range_handle = self.all_volumes
#         if vec_index.dtype == "bool":
#             vec = np.where(vec_index)[0]
#         else:
#             vec = vec_index.astype("uint")
#         handles = np.asarray(range_handle)[vec.astype("uint")].astype("uint")
#         return rng.Range(handles)
#
#     def set_data(self, name_tag, data, index_vec=np.array([]), range_el=None):
#         if range_el is None:
#             range_el = self.all_volumes
#         if index_vec.size > 0:
#             range_el = self.range_index(index_vec, range_el)
#         handle_tag = self.handleDic[name_tag]
#         self.mb.tag_set_data(handle_tag, range_el, data)
#
#     def init_tag(self, name_tag, dtype="int", entity_type=4):
#         # initialize a tag
#         # zeros nodes, edges, faces and volumes
#         # by default it zeros all geometric entities
#         if dtype == "int":
#             var_type = int
#         elif dtype == "float":
#             var_type = float
#         elif dtype == "bool":
#             var_type = bool
#         el = [[self.all_nodes], [self.all_edges], [self.all_faces], [self.all_volumes],
#               [self.all_nodes, self.all_edges, self.all_faces, self.all_volumes]]
#         range_temp = self.range_merge(*el[entity_type])
#         self.set_data(name_tag,data = np.zeros(len(range_temp)).astype(var_type),range_el = range_temp)
#
#     def check_range_by_dimm(self, handle):
#         # INPUT: handle or range
#         # OUTPUT: a vector with the same size as the input handle with the following classification
#         # 0 - nodes , 1 -edges, 2 - faces 3 - volumes, 4- meshset
#         handle_int = np.asarray(handle).astype("uint64")
#         type_list = np.array([self.mb.type_from_handle(el) for el in handle_int])
#         handle_classification = np.zeros(len(handle))
#         nodetype = type_list == types.MBVERTEX
#         edgetype = type_list == types.MBEDGE
#         facetype = (type_list == types.MBTRI) | (type_list == types.MBQUAD) | (type_list == types.MBPOLYGON)
#         volumetype = (type_list == types.MBTET) | (type_list == types.MBPYRAMID) | (type_list == types.MBPRISM) | \
#                       (type_list == types.MBKNIFE) | (type_list == types.MBHEX) | (type_list == types.MBPOLYHEDRON)
#         meshsettype = type_list == types.MBENTITYSET
#         handle_classification[nodetype] = 0
#         handle_classification[edgetype] = 1
#         handle_classification[facetype] = 2
#         handle_classification[volumetype] = 3
#         handle_classification[meshsettype] = 11
#         return handle_classification.astype("uint64")
#
#     def filter_range(self, handle, filter_vec):
#         # INPUT: handle or range
#         # OUTPUT: handles in the range of the dimension in args
#         # 0 - nodes , 1 -edges, 2 - faces 3 - volumes, 4- meshset
#         handle_int = np.asarray(handle).astype("uint64")
#         if filter_vec.dtype == bool:
#             return rng.Range(handle_int[filter_vec])
#         else:
#             return rng.Range(handle_int[filter_vec.astype("uint64")])
#
#     def filter_handle_by_dimension(self, handle, *args):
#         # INPUT: handle or range
#         # OUTPUT: handles in the range of the dimension in args
#         # 0 - nodes , 1 -edges, 2 - faces 3 - volumes, 4- meshset
#         handle_int = np.asarray(handle).astype("uint64")
#         vec_classification = self.check_range_by_dimm(handle)
#         test_elem = np.array([*args])
#         return rng.Range(handle_int[np.isin(vec_classification,test_elem)])
#
#     @staticmethod
#     def range_merge(*args):
#         range_merged = rng.Range()
#         for arg in args:
#                 range_merged.merge(arg)
#         return range_merged
#
#     def print(self, text=None):
#         m1 = self.mb.create_meshset()
#         self.mb.add_entities(m1, self.all_nodes)
#         m2 = self.mb.create_meshset()
#         self.mb.add_entities(m2, self.all_faces)
#         self.mb.remove_entities(m2, self.all_nodes)
#         m3 = self.mb.create_meshset()
#
#         self.mb.add_entities(m3, self.all_volumes)
#         m4 = self.mb.create_meshset()
#         self.mb.add_entities(m4, self.all_edges)
#         if text is None:
#             text = "output"
#         extension = ".vtk"
#         text1 = text + "-nodes" + extension
#         text2 = text + "-face" + extension
#         text3 = text + "-volume" + extension
#         text4 = text + "-edges" + extension
#         text5 = text + "-all" + extension
#         text6 = text + "-together" + extension
#         self.mb.write_file(text1, [m1])
#         self.mb.write_file(text2, [m2])
#         if self.dimension == 3:
#             self.mb.write_file(text3, [m3])
#
#         self.mb.write_file(text4, [m4])
#         self.mb.write_file(text5)
