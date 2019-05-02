"""
Module for implementation of multiscale mesh and CoarseVolumes objects functionalities
"""
#import time
import pdb
from . finescaleMesh import FineScaleMesh
from ..msCoarseningLib import algoritmo
from . meshComponents import MoabVariable
from . mscorePymoab import MsCoreMoab
from . meshComponentsMS import MoabVariableMS,  MeshEntitiesMS
from pymoab import core, types, rng
import numpy as np
import yaml


print('Initializing Finescale Mesh for Multiscale Methods')


class FineScaleMeshMS(FineScaleMesh):
    def __init__(self,mesh_file, dim = 3):
        super().__init__(mesh_file,dim)

        print("Creating Coarse Grid")
        self.coarse = MultiscaleCoarseGrid(self)
        # for i,el in zip(range(len(self.coarse_volumes)),self.coarse_volumes):
        #     el(i,self.general)

    def init_entities(self):
        self.nodes = MeshEntitiesMS(self.core, entity_type = "node")
        self.edges = MeshEntitiesMS(self.core, entity_type = "edges")
        self.faces = MeshEntitiesMS(self.core, entity_type = "faces")
        if self.dim == 3:
            self.volumes = MeshEntitiesMS(self.core, entity_type = "volumes")

    def init_variables(self):
        config = self.read_config('variable_settings.yml')

        nodes = config['nodes']
        edges = config['edges']
        faces = config['faces']
        volumes = config['volumes']
        not_empty = []
        parameters = [0,1]

        if nodes is not None:
            names = nodes.keys()
            for i in names:
                size = str(nodes[i]['data size'])
                format = nodes[i]['data format']
                command = 'self.' + i + ' = MoabVariableMS(self.core, data_size = ' + size + ', var_type = "nodes", data_format = ' + "'" + format + "'" + ', name_tag =' + "'" + i + "'" + ')'
                exec(command)
        if edges is not None:
            names = edges.keys()
            for i in names:
                size = str(edges[i]['data size'])
                format = edges[i]['data format']
                command = 'self.' + i + ' = MoabVariableMS(self.core, data_size = ' + size + ', var_type = "edges", data_format = ' + "'" + format + "'" + ', name_tag =' + "'" + i + "'" + ')'
                print(command)
                exec(command)
        if faces is not None:
            names = faces.keys()
            for i in names:
                size = str(faces[i]['data size'])
                format = faces[i]['data format']
                command = 'self.' + i + ' = MoabVariableMS(self.core, data_size = ' + size + ', var_type = "faces", data_format = ' + "'" + format + "'" + ', name_tag =' + "'" + i + "'" + ')'
                print(command)
                exec(command)
        if volumes is not None:
            names = volumes.keys()
            for i in names:
                size = str(volumes[i]['data size'])
                format = volumes[i]['data format']
                command = 'self.' + i + ' = MoabVariableMS(self.core, data_size = ' + size + ', var_type = "volumes", data_format = ' + "'" + format + "'" + ', name_tag =' + "'" + i + "'" + ')'
                print(command)
                exec(command)

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

    def init_partition_parallel(self):
        if self.dim == 3:
            partition = MoabVariable(self.core,data_size=1,var_type= "volumes",  data_format="int", name_tag="Parallel",
                                         data_density="sparse")

        elif self.dim == 2:
            partition = MoabVariable(self.core,data_size=1,var_type= "faces",  data_format="int", name_tag="Parallel", data_density="sparse")
        return partition

    def read_config(self, config_input):
        with open(config_input, 'r') as f:
            config_file = yaml.safe_load(f)
        return config_file


class CoarseVolume(FineScaleMeshMS):
    def __init__(self, father_core, dim, i, coarse_vec):
        self.dim = dim
        self.level = father_core.level + 1
        self.coarse_num = i

        print("Level {0} - Volume {1}".format(self.level,self.coarse_num))
        self.core = MsCoreMoab(father_core, i, coarse_vec)

        self.init_entities()
        self.init_variables()
        self.init_coarse_variables()
        self.macro_dim()

    def init_variables(self):
        #self.pressure = MoabVariableMS(self.core,data_size=1,var_type= "volumes",  data_format="float", name_tag="pressure", level=self.level, coarse_num=self.coarse_num)

        pass

    def __call__(self,i,general):
        self.nodes.enhance(i,general)
        self.edges.enhance(i,general)
        self.faces.enhance(i,general)
        if self.dim == 3:
            self.volumes.enhance(i,general)

    def init_coarse_variables(self):
        pass



class MultiscaleCoarseGrid(object):
    def __init__(self, M):
        print(M)
        self.mb = M.core.mb
        self.partition = M.init_partition()
        self.volumes = [CoarseVolume(M.core, M.dim, i, self.partition[:] == i) for i in range(self.partition[:].max()+1 )]
        self.num_coarse = len(self.volumes)
        self.num = {"nodes": 0, "node": 0, "edges": 1, "edge": 1, "faces": 2, "face": 2, "volumes": 3, "volume": 3,
                             0: 0, 1: 1, 2: 2, 3: 3}

        self.local_tag = [volume.core.handleDic[volume.core.id_name]  for volume in self.volumes]
        #self.local_tag =  [el.core.handleDic["LOCAL_ID_L" + str(el.core.level) + "-" + str(el.core.coarse_num)] for el in coarse_list]



        # for i,el in zip(range(len(self.coarse_volumes)),self.coarse_volumes):
        #     el(i,self.general)


        # for i,el in zip(range(len(self.coarse_volumes)),self.coarse_volumes):
        #     el(i,self.general)
        #self.num_coarse = len(coarse_list)
        self.find_coarse_neighbours()

        # self.local_tag = coarse_list[0].core.handleDic[coarse_list[0].core.id_name]
        self.global_tag = M.core.handleDic["GLOBAL_ID"]
        self.all_volumes = M.core.all_volumes
        self.all_faces = M.core.all_faces
        self.all_edges = M.core.all_edges
        self.all_nodes = M.core.all_nodes
        self.create_coarse_connectivities()

    def create_coarse_connectivities(self):
        self.connectivities = np.zeros((self.num_coarse,self.num_coarse))
        for x in range(self.num_coarse):
            for y in range(x+1,self.num_coarse):
                if not self.nodes_neighbors[x,y].empty():
                    self.connectivities[x,y]  =  True
                    self.connectivities[y,x]  =  True
    #
    def find_coarse_neighbours(self):
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
                self.nodes_neighbors[x,y] = rng.intersect(self.volumes[x].core.boundary_nodes, self.volumes[y].core.boundary_nodes)
                self.nodes_neighbors[y,x] = self.nodes_neighbors[x,y]
                temp = self.nodes_neighbors[x,y]
                [self.all_nodes_neighbors.insert(e) for e in temp]

                self.edges_neighbors[x,y] = rng.intersect(self.volumes[x].core.boundary_edges, self.volumes[y].core.boundary_edges)
                self.edges_neighbors[y,x] = self.edges_neighbors[x,y]
                temp = self.edges_neighbors[x,y]
                [self.all_edges_neighbors.insert(e) for e in temp]

                self.faces_neighbors[x,y] = rng.intersect(self.volumes[x].core.boundary_faces, self.volumes[y].core.boundary_faces)
                self.faces_neighbors[y,x] = self.faces_neighbors[x,y]
                temp = self.faces_neighbors[x,y]
                [self.all_faces_neighbors.insert(e) for e in temp]

                self.volumes_neighbors[x,y] = rng.intersect(self.volumes[x].core.boundary_volumes, self.volumes[y].core.boundary_volumes)
                self.volumes_neighbors[y,x] = self.volumes_neighbors[x,y]
                temp = self.volumes_neighbors[x,y]
                [self.all_volumes_neighbors.insert(e) for e in temp]
    #
    # def global_to_local_id(self,vec_range,element, target):
    #     flag = self.num[element]
    #     vec = self.create_range_vec(vec_range)
    #     if flag == 0:
    #         handle = self.range_index(vec, self.all_nodes)
    #     elif flag == 1:
    #         handle = self.range_index(vec, self.all_edges)
    #     elif flag == 2:
    #         handle = self.range_index(vec, self.all_faces)
    #     elif flag == 3:
    #         handle = self.range_index(vec, self.all_volumes)
    #     return self.mb.tag_get_data(self.local_tag[target],handle)
    #
    # def coarse_neighbours(self, x,y, element):
    #       # return self.read_data(self.global_tag, range_el = self.num[element])
    #       flag = self.num[element]
    #       if flag == 0:
    #           return self.mb.tag_get_data(self.global_tag, self.nodes_neighbors[x,y])
    #       elif flag == 1:
    #           return self.mb.tag_get_data(self.global_tag, self.edges_neighbors[x,y])
    #       elif flag == 2:
    #           return self.mb.tag_get_data(self.global_tag, self.faces_neighbors[x,y])
    #       elif flag == 3:
    #           return self.mb.tag_get_data(self.global_tag, self.volumes_neighbors[x,y])
    #
    # @property
    # def all_neighbors_nodes(self):
    #     return self.mb.tag_get_data(self.global_tag, self.all_nodes_neighbors)
    #
    # @property
    # def all_neighbors_edges(self):
    #     return self.mb.tag_get_data(self.global_tag, self.all_edges_neighbors)
    #
    # @property
    # def all_neighbors_faces(self):
    #     return self.mb.tag_get_data(self.global_tag, self.all_faces_neighbors)
    #
    # @property
    # def all_neighbors_volumes(self):
    #     return self.mb.tag_get_data(self.global_tag, self.all_volumes_neighbors)
    #
    # def create_range_vec(self, index):
    #     range_vec = None
    #     if isinstance(index, int) or isinstance(index, np.integer):
    #         range_vec = np.array([index]).astype("uint")
    #     elif isinstance(index, np.ndarray):
    #         if index.dtype == "bool":
    #             range_vec = np.where(index)[0]
    #         else:
    #             range_vec = index
    #     elif isinstance(index, slice):
    #         start = index.start
    #         stop = index.stop
    #         step = index.step
    #         if start is None:
    #             start = 0
    #         if stop is None:
    #             stop = len(self)
    #         if step is None:
    #             step = 1
    #         if start < 0:
    #             start = len(self) + start + 1
    #         if stop < 0:
    #             stop = len(self) + stop + 1
    #         range_vec = np.arange(start, stop, step).astype('uint')
    #     elif isinstance(index, list):
    #         range_vec = np.array(index)
    #     return range_vec
    #
    # def read_data(self, name_tag, index_vec = np.array([]), range_el = None):
    #     if range_el is None:
    #         range_el = self.all_volumes
    #     if index_vec.size > 0:
    #         range_el = self.range_index(index_vec,range_el)
    #     try:
    #         handle_tag = self.handleDic[name_tag]
    #         return self.mb.tag_get_data(handle_tag, range_el)
    #     except KeyError:
    #         print("Tag not found")
    #
    # def range_index(self, vec_index, range_handle = None):
    #     if range_handle is None:
    #         range_handle = self.all_volumes
    #     if vec_index.dtype == "bool":
    #         vec = np.where(vec_index)[0]
    #     else:
    #         vec = vec_index.astype("uint")
    #     handles = np.asarray(range_handle)[vec.astype("uint")].astype("uint")
    #     return rng.Range(handles)
