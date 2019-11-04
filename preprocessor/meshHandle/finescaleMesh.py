"""
Module for management of fine scale mesh
"""
import time
import pdb
import numpy as np
from ..meshHandle.meshComponents import MoabVariable
from ..meshHandle.configTools.configClass import variableInit
from math import sqrt
from pymoab import core, types, rng, topo_util
from . corePymoab import CoreMoab
from . meshComponents import MeshEntities, MoabVariable
import yaml

print('Standard fine-scale mesh loaded: No multiscale components available')


class FineScaleMesh:
    def __init__(self, mesh_file, dim=3, var_config=None):
        self.var_config = var_config
        self.dim = dim
        self.core = CoreMoab(mesh_file, dim)
        if mesh_file is not None:
            self.run()

    def run(self):
        self.init_entities()
        self.init_variables()
        self.macro_dim()
        self.init_dimmension()

    def init_variables(self):
        if self.var_config is None:
            self.var_config = variableInit()
        for command in self.var_config.get_var(self.core.level):
            exec(command)

    def init_entities(self):
        self.nodes = MeshEntities(self.core, entity_type = "nodes")
        self.edges = MeshEntities(self.core, entity_type="edges")
        self.faces = MeshEntities(self.core, entity_type = "faces")
        if self.dim == 3:
            self.volumes = MeshEntities(self.core, entity_type = "volumes")

    def __len__(self):
        if self.dim == 3:
            return len(self.volumes)
        elif self.dim ==2:
            return len(self.faces)

    def macro_dim(self):
        # coords = self.core.mb.get_coords(self.core.all_nodes).reshape(len(self.core.all_nodes),3)
        min_coord = self.nodes.coords[:].min(axis = 0)
        max_coord  = self.nodes.coords[:].max(axis = 0)
        self.rx = (min_coord[0], max_coord[0])
        self.ry = (min_coord[1], max_coord[1])
        self.rz= (min_coord[2], max_coord[2])

    def init_dimmension(self):
        # center = MoabVariable(self.core, data_size=3, var_type=0, data_format='float', name_tag='CENTER', data_density='dense')
        # import pdb; pdb.set_trace()
        pass

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


    def init_volume(self):
        mat = np.zeros(self.core.all_volumes.size())
        index = 0
        for vol in  self.core.all_volumes:
            mat[index] = self.get_volume(vol)
            index += 1
        self.core.create_tag_handle("VOLUME",1)
        #pdb.set_trace()
        self.core.set_data("VOLUME",mat)
        #volTag = self.mb.tag_get_handle("VOLUME", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
        #self.mb.tag_set_data(volTag, self.all_volumes, mat)
        #return volTag

    def get_centroid(self, entity):
        verts = self.core.mb.get_connectivity(entity)
        coords = np.array([self.core.mb.get_coords([vert]) for vert in verts])
        qtd_pts = len(verts)
        #print qtd_pts, 'qtd_pts'
        coords = np.reshape(coords, (qtd_pts, 3))
        pseudo_cent = sum(coords)/qtd_pts
        return pseudo_cent

    def get_volume(self,entity):
        #input: entity tag
        #ouput: volume of a entity
        verts = self.core.mb.get_connectivity(entity)
        coords = np.array([self.core.mb.get_coords([vert]) for vert in verts])
        qtd_pts = len(verts)
        if qtd_pts == 4:
            pass
            vect_1 = coords[1] - coords[0]
            vect_2 = coords[2] - coords[0]
            vect_3 = coords[3] - coords[0]
            vol_eval = abs(np.dot(np.cross(vect_1, vect_2), vect_3)) / 6.0
        elif qtd_pts == 8:
            pass
            #SEGUNDA ATIVIDADE PRA RENATINHA
            #CALCULAR O VOLUME DO HEXAEDRO IERREGULAR DADO OS 8 PONTOS

            # pass
            # #pdb.set_trace()
            # vect_1 = coords[7] - coords[0]
            # vect_2 = coords[1] - coords[0]
            # vect_3 = coords[3] - coords[5]
            # vect_4 = coords[4] - coords[6]
            # vect_5 = coords[5] - coords[0]
            # vect_6 = coords[2] - coords[0]
            # vect_7 = coords[6] - coords[3]
            # D1 = np.linalg.det(np.array([vect_1, vect_2, vect_3]))
            # D2 = np.linalg.det(np.array([vect_1, vect_4, vect_5]))
            # D3 = np.linalg.det(np.array([vect_1, vect_6, vect_7]))
            # vol_eval = ((D1+D2+D3)/2)
            vol_eval = 1
            return vol_eval
        else:
            vol_eval  = 0
        return vol_eval

    def get_tetra_volume(self, tet_nodes):
        vect_1 = tet_nodes[1] - tet_nodes[0]
        vect_2 = tet_nodes[2] - tet_nodes[0]
        vect_3 = tet_nodes[3] - tet_nodes[0]
        vol_eval = abs(np.dot(np.cross(vect_1, vect_2), vect_3))/1
        return vol_eval

    def get_piram_volume(self,tet_nodes):
        pass

    @staticmethod
    def point_distance(coords_1, coords_2):
        dist_vector = coords_1 - coords_2
        distance = sqrt(np.dot(dist_vector, dist_vector))
        return distance
