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
        # import pdb; pdb.set_trace()
        self.coarse = MultiscaleCoarseGrid(self)
        self.enhance_entities()


    def enhance_entities(self):
        for i,el in zip(range(len(self.coarse.elements)),self.coarse.elements):
            el(i,self.coarse)

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
        self.nodes.enhance(i, general)
        self.edges.enhance(i, general)
        self.faces.enhance(i, general)
        if self.dim == 3:
            self.volumes.enhance(i,general)

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
        tmp = self.dic[item]
        el_list = rng.Range()
        if not isinstance(item, int):
            for e in tmp:
                el_list.insert(e[0])
            return self.fun(self.tag, el_list)
        else:
            return self.fun(self.tag, tmp)



class MultiscaleCoarseGrid(object):
    def __init__(self, M):
        self.mb = M.core.mb
        self.partition = M.init_partition()
        self.elements = [CoarseVolume(M.core, M.dim, i, self.partition[:] == i) for i in range(self.partition[:].max()+1 )]
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



    def find_coarse_neighbours(self):
        self.connectivities = np.zeros((self.num_coarse,self.num_coarse+1 ,3)).astype('bool')
        self.nodes_neighbors  = np.zeros((self.num_coarse,self.num_coarse+1), dtype = object)
        self.edges_neighbors  = np.zeros((self.num_coarse,self.num_coarse+1), dtype = object)
        self.faces_neighbors  = np.zeros((self.num_coarse,self.num_coarse+1), dtype = object)
        self.nodes_neighbors[:], self.edges_neighbors[:], self.faces_neighbors[:] = None , None, None
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
                    self.nodes_neighbors[x,y], self.nodes_neighbors[y,x],= node_count ,node_count
                    self.connectivities[x, y, 0],self.connectivities[y, x, 0] = True, True
                    node_count += 1
                    [self.all_nodes_neighbors.insert(e) for e in node_intersect]
                edges_intersect = rng.intersect(self.elements[x].core.boundary_edges, self.elements[y].core.boundary_edges)
                if not edges_intersect.empty():
                    self._edges.append(edges_intersect)
                    # self._edges = np.append(self._edges,edges_intersect)
                    self.edges_neighbors[x,y], self.edges_neighbors[y,x]= edge_count ,edge_count
                    self.connectivities[x, y, 1], self.connectivities[y, x, 1] =  True, True
                    edge_count += 1
                    [self.all_edges_neighbors.insert(e) for e in edges_intersect]
                faces_intersect = rng.intersect(self.elements[x].core.boundary_faces, self.elements[y].core.boundary_faces)
                if not faces_intersect.empty():
                    self._faces.append(faces_intersect)
                    #self._faces = np.append(self._faces,faces_intersect)
                    self.faces_neighbors[x,y], self.faces_neighbors[y,x]= face_count ,face_count
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
                self.nodes_neighbors[x, -1] = node_count
                self.connectivities[x, -1, 0] = True
                node_count += 1
            edge_intersect = rng.subtract(self.elements[x].core.boundary_edges, self.all_edges_neighbors)
            if not edge_intersect.empty():
                self._edges.append(edge_intersect)
                self.edges_neighbors[x, -1] = edge_count
                self.connectivities[x, -1, 1] = True
                edge_count += 1
            face_intersect = rng.subtract(self.elements[x].core.boundary_faces, self.all_faces_neighbors)
            if not face_intersect.empty():
                self._faces.append(face_intersect)
                self.faces_neighbors[x, -1] = face_count
                self.connectivities[x, -1, 2] = True
                face_count += 1

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
