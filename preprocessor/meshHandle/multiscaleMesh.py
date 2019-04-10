"""
Module for implementation of multiscale mesh and CoarseVolumes objects functionalities
"""
#import time
import pdb
from . finescaleMesh import FineScaleMesh
from ..msCoarseningLib import algoritmo
from . meshComponents import MoabVariable
from . mscorePymoab import MsCoreMoab
from . meshComponentsMS import MultiscaleMeshEntities, MoabVariableMS,  MeshEntitiesMS
from pymoab import core, types, rng
import yaml


print('Initializing Finescale Mesh for Multiscale Methods')


class FineScaleMeshMS(FineScaleMesh):
    def __init__(self,mesh_file, dim = 3):
        super().__init__(mesh_file,dim)
        self.partition = self.init_partition()
        self.coarse_volumes = [CoarseVolume(self.core, self.dim, i, self.partition[:] == i) for i in range(self.partition[:].max()+1 )]
        print("Creating object general")
        self.general = MultiscaleMeshEntities(self.core,self.coarse_volumes)
        for i,el in zip(range(len(self.coarse_volumes)),self.coarse_volumes):
            el(i,self.general)

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
                exec(command)
        if faces is not None:
            names = faces.keys()
            for i in names:
                size = str(faces[i]['data size'])
                format = faces[i]['data format']
                command = 'self.' + i + ' = MoabVariableMS(self.core, data_size = ' + size + ', var_type = "faces", data_format = ' + "'" + format + "'" + ', name_tag =' + "'" + i + "'" + ')'
                exec(command)
        if volumes is not None:
            names = volumes.keys()
            for i in names:
                size = str(volumes[i]['data size'])
                format = volumes[i]['data format']
                command = 'self.' + i + ' = MoabVariableMS(self.core, data_size = ' + size + ', var_type = "volumes", data_format = ' + "'" + format + "'" + ', name_tag =' + "'" + i + "'" + ')'
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

            # partition[:]
            # [partition[:],coarse_center]  = getattr(msCoarseningLib.algoritmo, name_function)(self.volumes.center[:],
            #            len(self), self.rx, self.ry, self.rz,*used_attributes)
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
        pass

    def __call__(self,i,general):
        self.nodes.enhance(i,general)
        self.edges.enhance(i,general)
        self.faces.enhance(i,general)
        if self.dim == 3:
            self.volumes.enhance(i,general)

    def init_coarse_variables(self):
        pass
