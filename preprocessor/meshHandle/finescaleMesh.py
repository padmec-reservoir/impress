"""
Module for management of fine scale mesh
"""
import numpy as np
from ..meshHandle.serialization import IMPRESSPickler, IMPRESSUnpickler
from ..meshHandle.meshComponents import MoabVariable
from ..meshHandle.configTools.configClass import variableInit
from math import sqrt
from pymoab import core, types, rng, topo_util
from . corePymoab import CoreMoab
from . meshComponents import MeshEntities, MoabVariable
import yaml
import pickle

print('Standard fine-scale mesh loaded: No multiscale components available')

class FineScaleMesh:
    def __init__(self, mesh_file, dim=3, var_config=None, load = False):
        self.load = load
        self.mesh_file = mesh_file
        self.var_config = var_config
        self.dim = dim
        self.core = CoreMoab(mesh_file, dim)
        if mesh_file is not None:
            self.run()

    def run(self):
        self.init_entities()
        if not self.load:
            self.init_variables()
        else:
            self.load_variables()
        self.macro_dim()

    def create_variable(self, name_tag, var_type="volumes", data_size=1, data_format="float", data_density="sparse",
                 entity_index=None, create = True):
        var = MoabVariable(self.core, data_size= data_size, var_type = var_type, data_format = data_format, name_tag = name_tag, data_density = data_density, entity_index = entity_index, create = create)
        exec(f'self.{name_tag} = var')
        self.var_handle_list.append(var)
        return var

    def save_variables(self, name_file):
        self.core.mb.write_file('saves/'+name_file+'.h5m')
        file = open('saves/'+name_file+'.imp', 'wb')
        pickle.dump([(tags.name_tag, tags.var_type, tags.data_size, tags.data_format, tags.data_density) for tags in self.var_handle_list], file)
        file.close()

    def dump(self, file_name):
        file = open(file_name, "wb+")
        pic = IMPRESSPickler(file)
        pic.dump(self)

    def init_variables(self):
        self.var_handle_list = []
        if self.var_config is None:
            self.var_config = variableInit()
        for command in self.var_config.get_var(self.core.level):
            exec(command)

    def load_variables(self):
        self.var_handle_list = []
        file = open(self.mesh_file.split('.')[0]+'.imp', 'rb')
        tag_list = pickle.load(file)
        file.close()
        for tags in tag_list:
            self.create_variable(name_tag = tags[0], var_type = tags[1], data_size = tags[2], data_format = tags[3], data_density = tags[4], create = False)

    def init_entities(self):
        self.nodes = MeshEntities(self.core, entity_type = "nodes")
        self.edges = MeshEntities(self.core, entity_type="edges")
        self.faces = MeshEntities(self.core, entity_type = "faces")
        if self.dim == 3:
            self.volumes = MeshEntities(self.core, entity_type = "volumes")

    def __len__(self):
        if self.dim == 3:
            return len(self.volumes)
        elif self.dim == 2:
            return len(self.faces)

    def load_array(self, type = None, array = None):
        if type == None:
            self.nodes.load_array(array)
            self.edges.load_array(array)
            self.faces.load_array(array)
            self.volumes.load_array(array)
        if type == 'nodes' or type == 0:
            self.nodes.load_array(array)
        if type == 'edges' or type == 1:
            self.edges.load_array(array)
        if type == 'faces' or type == 2:
            self.faces.load_array(array)
        if type == 'volumes' or type == 3:
            self.volumes.load_array(array)

    def to_moab(self):
        for vars in self.var_handle_list:
            vars.to_moab()

    def to_numpy(self):
        for variables in self.var_handle_list:
            variables.to_numpy()

    def unload_array(self, type = None, array = None):
        if type == None:
            self.nodes.unload_array(array)
            self.edges.unload_array(array)
            self.faces.unload_array(array)
            self.volumes.unload_array(array)
        if type == 'nodes' or type == 0:
            self.nodes.unload_array(array)
        if type == 'edges' or type == 1:
            self.edges.unload_array(array)
        if type == 'faces' or type == 2:
            self.faces.unload_array(array)
        if type == 'volumes' or type == 3:
            self.volumes.unload_array(array)

    def macro_dim(self):
        min_coord = self.nodes.coords[:].min(axis = 0)
        max_coord  = self.nodes.coords[:].max(axis = 0)
        self.rx = (min_coord[0], max_coord[0])
        self.ry = (min_coord[1], max_coord[1])
        self.rz= (min_coord[2], max_coord[2])

    def init_center(self):
        self.core.create_tag_handle('CENTER',3)
        #centro dos volumes
        centers = np.zeros((len(self.core.all_volumes),3)).astype('float')
        index = 0
        for volume in self.core.all_volumes:
            centers[index] = self.get_centroid(volume)
            index += 1
        self.core.set_data("CENTER",centers)
        #centro das faces
        centers = np.zeros((len(self.core.all_faces), 3)).astype('float')
        index = 0
        for face in self.core.all_faces:
            centers[index] = self.get_centroid(face)
            index += 1
        self.core.set_data("CENTER", centers, range_el = self.core.all_faces)

        #centro das arestas
        centers = np.zeros((len(self.core.all_edges), 3)).astype('float')
        index = 0
        for edge in self.core.all_edges:
            centers[index] = self.get_centroid(edge)
            index += 1
        self.core.set_data("CENTER", centers, range_el = self.core.all_edges)

    def init_normal(self):
        self.core.create_tag_handle('NORMAL', 3)
        normal = np.zeros((len(self.core.all_faces), 3)).astype('float')
        index = 0
        for face in self.core.all_faces:
            verts = self.core.mb.get_connectivity(face)
            coords = np.array([self.core.mb.get_coords([vert]) for vert in verts])
            vec1 = coords[1] - coords[0]
            vec2 = coords[2] - coords[0]
            cross = np.cross(vec1,vec2)
            normal[index] = cross/np.linalg.norm(cross)
            index += 1
        self.core.set_data("NORMAL", normal, range_el=self.core.all_faces)

    def get_centroid(self, entity):
        verts = self.core.mb.get_connectivity(entity)
        coords = np.array([self.core.mb.get_coords([vert]) for vert in verts])
        qtd_pts = len(verts)
        coords = np.reshape(coords, (qtd_pts, 3))
        pseudo_cent = sum(coords)/qtd_pts
        return pseudo_cent
