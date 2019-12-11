"""
Use of Pymoab methods to read and manage the input mesh
"""
from pymoab import core, types, rng, topo_util
from pymoab import skinner as sk
import numpy as np
import yaml
import time



class CoreMoab:
    def __init__(self, mesh_file=None, dim=3):
        self.dimension = dim
        self.mb = core.Core()
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        if mesh_file is not None:
            self.mb.load_file(mesh_file)
            self.run()

    def run(self):
        self.root_set = self.mb.get_root_set()
        self.father_root_set = self.root_set
        self.level = 0
        self.all_volumes = self.mb.get_entities_by_dimension(0, 3)
        self.all_nodes = self.mb.get_entities_by_dimension(0, 0)
        self.mtu.construct_aentities(self.all_nodes)
        #self.mtu.construct_aentities(self.all_volumes)
        self.all_faces = self.mb.get_entities_by_dimension(0, 2)
        self.all_edges = self.mb.get_entities_by_dimension(0, 1)
        self.handleDic = {}
        [self.boundary_nodes, self.boundary_edges, self.boundary_faces, self.boundary_volumes] = self.skinner_operation()
        self.internal_nodes = rng.subtract(self.all_nodes, self.boundary_nodes)
        self.internal_edges = rng.subtract(self.all_edges, self.boundary_edges)
        self.internal_faces = rng.subtract(self.all_faces, self.boundary_faces)
        self.internal_volumes = rng.subtract(self.all_volumes, self.boundary_volumes)
        self.id_name = "GLOBAL_ID"
        self.father_id_name = self.id_name
        self.init_id()
        #self.check_integrity()
        # self.create_id_visualization()
        # self.flag_dic = {}
        [self.flag_list, self.flag_dic] = self.read_flags()

        # self.create_flag_visualization()

        # swtich on/off
        self.parallel_meshset = self.create_parallel_meshset()
        self.create_parallel_tag()

    def init_id(self):
        self.create_tag_handle(self.id_name, data_size = 1, data_text = "int", data_density = "sparse")
        # create volume ids
        self.set_data(self.id_name, np.arange(len(self.all_volumes)))
        # create face ids
        self.set_data(self.id_name, np.arange(len(self.all_faces)), range_el=self.all_faces)
        # create edges ids
        self.set_data(self.id_name, np.arange(len(self.all_edges)), range_el=self.all_edges)
        # create nodes ids
        self.set_data(self.id_name, np.arange(len(self.all_nodes)), range_el=self.all_nodes)

    def skinner_operation(self):
        self.skin = sk.Skinner(self.mb)
        print("Entering skinner test")

        if self.dimension == 3:
            faces_on_skin_handles = self.skin.find_skin(self.root_set, self.all_volumes)
            edges_on_skin_handles = self.mtu.get_bridge_adjacencies(faces_on_skin_handles, 2, 1)
            nodes_on_skin_handles = self.mtu.get_bridge_adjacencies(faces_on_skin_handles, 2, 0)
            volumes_on_skin_handles = self.mtu.get_bridge_adjacencies(faces_on_skin_handles, 0, 3)
        elif self.dimension == 2:
            edges_on_skin_handles = self.skin.find_skin(self.root_set, self.all_faces)
            nodes_on_skin_handles = self.mtu.get_bridge_adjacencies(edges_on_skin_handles, 1, 0)
            faces_on_skin_handles = self.mtu.get_bridge_adjacencies(edges_on_skin_handles, 0, 2)
            volumes_on_skin_handles = self.mtu.get_bridge_adjacencies(edges_on_skin_handles, 0, 3) #empty

        print("Skinning Operation Successful")
        return [nodes_on_skin_handles, edges_on_skin_handles, faces_on_skin_handles, volumes_on_skin_handles]


    def check_integrity(self):
        # check if the mesh contains
        check_list = [len(self.all_nodes), len(self.all_edges), len(self.all_faces), len(self.all_volumes)]
        list_words = ['Nodes', "Edges", "Faces", "Volumes"]
        print("Checking mesh integrity:")
        index = 0
        if self.dimension == 2:
            list_words = list_words[:-1]
            check_list = check_list[:-1]
        for entity in check_list:
            if entity > 0:
                print(list_words[index] + " successfully imported")
            else:
                print("------------------------------\nError creating \n" +
                      list_words[index] + " was not imported")
            index += 1

    def create_parallel_meshset(self):
        partition_volumes = []
        parallel_tag = []
        try:
            parallel_tag = self.mb.tag_get_handle("PARALLEL_PARTITION")
            flag = False
        except:
            print("Parallel Partition Tag not found \nAborting creating parallel partition entities.")
            flag = True
        if not flag:
            print("Parallel Partition Tag detected \nCreating parallel partitions entities")
            self.handleDic["PARALLEL_PARTITION"] = parallel_tag
            parallel_sets = self.mb.get_entities_by_type_and_tag(
                0, types.MBENTITYSET, np.array(
                (parallel_tag,)), np.array((None,)))
            self.create_tag_handle("Parallel", data_size=1, data_text="int")
            # partition_volumes = []
            for set_el in parallel_sets:
                num_tag = self.read_data("PARALLEL_PARTITION", range_el=set_el)[0, 0]
                # self.set_data("Parallel", np.ones(len(set_el))*num_tag,  range_el=set_el)
                list_entity = [self.mb.get_entities_by_dimension(set_el, 0), self.mb.get_entities_by_dimension(set_el, 1),
                               self.mb.get_entities_by_dimension(set_el, 2), self.mb.get_entities_by_dimension(set_el, 3)]
                # print([num_tag, list_entity])
                partition_volumes.append(list_entity)
        return partition_volumes

    def create_parallel_tag(self):
        k = 0
        for sets in self.parallel_meshset:
            for dim in sets:
                if len(dim) != 0:
                    self.set_data("Parallel", k * np.ones(len(dim)).astype(int), range_el=dim)
            k += 1

    def read_flags(self):
        physical_tag = self.mb.tag_get_handle("MATERIAL_SET")
        physical_sets = self.mb.get_entities_by_type_and_tag(
            0, types.MBENTITYSET, np.array(
            (physical_tag,)), np.array((None,)))
        self.handleDic["MATERIAL_SET"] = physical_tag
        flag_list = np.array([])
        flag_dic = {}
        for set in physical_sets:
            bc_flag = self.read_data("MATERIAL_SET", range_el=set)[0, 0]
            flag_list = np.append(flag_list, bc_flag)
            list_entity = [self.mb.get_entities_by_dimension(set, 0), self.mb.get_entities_by_dimension(set, 1),
                           self.mb.get_entities_by_dimension(set, 2), self.mb.get_entities_by_dimension(set, 3)]
            flag_dic[bc_flag] = list_entity
        return np.sort(flag_list), flag_dic

    # def access_meshset(self, handle):
    #     # returns the entities contained inside a give meshset handle
    #     # ie: for a meshset handle the entities inside are returned
    #     temp_range = []
    #     for el in range(self.dimension+1):
    #         sub_el = (self.mb.get_entities_by_dimension(handle, el))
    #         temp_range.append(sub_el)
    #     temp_range.append(self.mb.get_entities_by_dimension(handle, 11))
    #     return temp_range

    # def access_handle(self,handle):
    #     # input: range of handles of different dimensions
    #     # returns all entities with d-1 dimension the comprises the given range
    #     # ie: for a volume, the faces, for a face the edges and for an edge the points.
    #     #
    #     vecdim = self.check_range_by_dimm(handle)
    #     # pdb.set_trace()
    #     all_adj = np.array([np.array(self.mb.get_adjacencies(el_handle, dim-1)) for dim, el_handle in zip(vecdim,handle)])
    #     #unique_adj = np.unique(np.ma.concatenate(all_adj)).astype("uint64")
    #     # print(handle,all_adj)
    #     unique_adj = np.unique(np.concatenate(all_adj)).astype("uint64")
    #     return rng.Range(unique_adj)

    def create_tag_handle(self, name_tag, data_size, data_text="float", data_density="dense"):
        if data_density == "dense":
            data_density = types.MB_TAG_DENSE
        elif data_density == "sparse":
            data_density = types.MB_TAG_SPARSE
        elif data_density == "bit":
            data_density = types.MB_TAG_BIT
        else:
            print("Please define a valid tag type")
        if data_text == 'float':
            data_type = types.MB_TYPE_DOUBLE
        elif data_text == "int":
            data_type = types.MB_TYPE_INTEGER
        elif data_text == "bool":
            data_type = types.MB_TYPE_BIT
        try:
            handle = self.handleDic[name_tag]
        except KeyError:
            handle = self.mb.tag_get_handle(name_tag, data_size, data_type, data_density, True)
            self.handleDic[name_tag] = handle

    def read_data(self, name_tag, index_vec = np.array([]), range_el = None):
        if range_el is None:
            range_el = self.all_volumes
        if index_vec.size > 0:
            range_el = range_el[index_vec]
        try:
            handle_tag = self.handleDic[name_tag]
            return self.mb.tag_get_data(handle_tag, range_el)
        except KeyError:
            print("Tag not found")

    # def range_index(self, vec_index, range_handle=None):
    #     if range_handle is None:
    #         range_handle = self.all_volumes
    #     if vec_index.dtype == "bool":
    #         vec = np.where(vec_index)[0]
    #     else:
    #         vec = vec_index.astype("uint")
    #     handles = np.asarray(range_handle)[vec.astype("uint")].astype("uint")
    #     return rng.Range(handles)

    def set_data(self, name_tag, data, index_vec=np.array([]), range_el=None):
        if range_el is None:
            range_el = self.all_volumes
        if index_vec.size > 0:
            range_el = range_el[index_vec]
        handle_tag = self.handleDic[name_tag]
        self.mb.tag_set_data(handle_tag, range_el.get_array(), data)

    def check_range_by_dimm(self, handle):
        # INPUT: handle or range
        # OUTPUT: a vector with the same size as the input handle with the following classification
        # 0 - nodes , 1 -edges, 2 - faces 3 - volumes, 4- meshset
        handle_int = np.asarray(handle).astype("uint64")
        type_list = np.array([self.mb.type_from_handle(el) for el in handle_int])
        handle_classification = np.zeros(len(handle))
        nodetype = type_list == types.MBVERTEX
        edgetype = type_list == types.MBEDGE
        facetype = (type_list == types.MBTRI) | (type_list == types.MBQUAD) | (type_list == types.MBPOLYGON)
        volumetype = (type_list == types.MBTET) | (type_list == types.MBPYRAMID) | (type_list == types.MBPRISM) | \
                      (type_list == types.MBKNIFE) | (type_list == types.MBHEX) | (type_list == types.MBPOLYHEDRON)
        meshsettype = type_list == types.MBENTITYSET
        handle_classification[nodetype] = 0
        handle_classification[edgetype] = 1
        handle_classification[facetype] = 2
        handle_classification[volumetype] = 3
        handle_classification[meshsettype] = 11
        return handle_classification.astype("uint64")

    def filter_range(self, handle, filter_vec):
        # INPUT: handle or range
        # OUTPUT: handles in the range of the dimension in args
        # 0 - nodes , 1 -edges, 2 - faces 3 - volumes, 4- meshset
        handle_int = np.asarray(handle).astype("uint64")
        if filter_vec.dtype == bool:
            return rng.Range(handle_int[filter_vec])
        else:
            return rng.Range(handle_int[filter_vec.astype("uint64")])

    def filter_handle_by_dimension(self, handle, *args):
        # INPUT: handle or range
        # OUTPUT: handles in the range of the dimension in args
        # 0 - nodes , 1 -edges, 2 - faces 3 - volumes, 4- meshset
        handle_int = np.asarray(handle).astype("uint64")
        vec_classification = self.check_range_by_dimm(handle)
        test_elem = np.array([*args])
        return rng.Range(handle_int[np.isin(vec_classification,test_elem)])

    @staticmethod
    def range_merge(*args):
        range_merged = rng.Range()
        for arg in args:
                range_merged.merge(arg)
        return range_merged

    def print(self, file=None, extension=".h5m", case = None,  config_input="input_cards/print_settings.yml"):
        if case is None:
            case = ''
        text =  file
        folder = "results/" + case
        with open(config_input, 'r') as f:
            data = yaml.safe_load(f)
        nodes = data['nodes']
        edges = data['edges']
        faces = data['faces']
        volumes = data['volumes']
        all_entities = data['all entities']
        m1 = self.mb.create_meshset()
        self.mb.add_entities(m1, self.all_nodes)
        m2 = self.mb.create_meshset()
        self.mb.add_entities(m2, self.all_faces)
        self.mb.remove_entities(m2, self.all_nodes)
        m3 = self.mb.create_meshset()
        self.mb.add_entities(m3, self.all_volumes)
        m4 = self.mb.create_meshset()
        self.mb.add_entities(m4, self.all_edges)
        if text is None:
            text = "output"
        text1 = text + "-nodes" + extension
        text2 = text + "-face" + extension
        text3 = text + "-volume" + extension
        text4 = text + "-edges" + extension
        text5 = text + "-all" + extension
        if folder is not None:
            import os
            if not os.path.exists(folder):
                os.mkdir(folder)
            text = folder + '/' + text
            text1 = folder + '/' + text1
            text2 = folder + '/' + text2
            text3 = folder + '/' + text3
            text4 = folder + '/' + text4
            text5 = folder + '/' + text5
        if nodes != 0:
            self.mb.write_file(text1, [m1])
        if faces != 0:
            self.mb.write_file(text2, [m2])
        if self.dimension == 3 and volumes != 0:
            self.mb.write_file(text3, [m3])
        if edges != 0:
            self.mb.write_file(text4, [m4])
        if all_entities != 0:
            self.mb.write_file(text5)
