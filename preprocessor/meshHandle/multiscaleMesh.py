"""
Module for implementation of multiscale mesh and CoarseVolumes objects functionalities
"""
#import time
#import pdb
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
        #self.alma = MoabVariableMS(self.core,data_size=1,var_type= "volumes",  data_format="int", name_tag="alma")
        #self.ama = MoabVariableMS(self.core,data_size=1,var_type= "faces",  data_format="float", name_tag="ama",data_density="sparse")
        #self.arma = MoabVariableMS(self.core,data_size=3,var_type= "edges",  data_format="float", name_tag="arma", data_density="sparse")
        self.permeability = MoabVariableMS(self.core,data_size=1,var_type= "volumes",  data_format="int", name_tag="permeability")
        self.pressure = MoabVariableMS(self.core,data_size=1,var_type= "volumes",  data_format="float", name_tag="pressure")
        self.erro = MoabVariableMS(self.core,data_size=1,var_type= "volumes",  data_format="float", name_tag="erro")


    def init_partition(self):
        config = self.read_config()
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
                key = "Coarsening_" + particionador_type + "_Input"
                specific_attributes = config.items(key)
                used_attributes = []
                for at in specific_attributes:
                    used_attributes.append(float(at[1]))
                [partition[:],coarse_center]  = getattr(msCoarseningLib.algoritmo, name_function)(self.faces.center[:],
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

    def read_config(self, config_input="msCoarse.yml"):
        with open("msCoarse.yml", 'r') as f:
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
        #self.lama = MoabVariableMS(self.core,data_size=1,var_type= "faces",  data_format="int", name_tag="lama", level=self.level, coarse_num=self.coarse_num)
        self.pressure_coarse = MoabVariableMS(self.core,data_size=1,var_type= "volumes",  data_format="float", name_tag="pressure_coarse", level=self.level, coarse_num=self.coarse_num)
        self.permeability_coarse = MoabVariableMS(self.core,data_size=1,var_type= "volumes",  data_format="float", name_tag="permeability_coarse", level=self.level, coarse_num=self.coarse_num)
