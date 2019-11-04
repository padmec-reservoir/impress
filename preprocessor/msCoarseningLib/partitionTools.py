

import numpy as np
import yaml
import os
from ..meshHandle.finescaleMesh import FineScaleMesh
from ..meshHandle.configTools.configClass import variableInit
from pymoab import types
import copy


class partitionManager(object):
    def __init__(self, M, config_object):
        scheme = config_object.tree['Scheme']
        if scheme == 'smart':
            self.func = self.smart
            if config_object.tree['Smart']['path'] == 'default':
                file_path = os.getcwd() + '/mesh/coarse/' + \
                    config_object.tree['Smart']['file']
            else:
                file_path = config_object.tree['Smart']['path'] + '/' + \
                    config_object.tree['Smart']['file']
            self.arg = [file_path]
            self.partitioner = smartPartition(M)
        elif scheme == 'simple':
            self.func = self.simple
            self.partitioner = simplePartition(M)
            self.arg = config_object.tree['Simple']['nx'], config_object.tree['Simple']['ny'], config_object.tree['Simple']['nz']
        elif scheme == 'parallel':
            print('Looking for Parallel Partition')
            self.func = self.parallel()
        else:
            print('Scheme ' + scheme + ' not defined.')

    def __call__(self):
        return self.run()

    def run(self, *args):
        if len(args) == 0:
            return self.func(*self.arg)
        else:
            return self.func(*args)

    def parallel(self):
        return ['parallel', None]

    def smart(self, file_path):
        return self.partitioner(file_path)

    def simple(self, nx, ny, nz):
        return self.partitioner(nx, ny, nz)


class smartPartition(object):
    def __init__(self,M):
        self.M = M

    def __call__(self, file_path):
        return self.run(file_path)

    def run(self, file_path):
        print('Creating Coarse Scale Forming Primal Grid')
        self.variable_entries = variableInit(empty=True)
        self.variable_entries.add_entry('pressure', 'volumes', 1, 'float')
        self.primal = FineScaleMesh(file_path, var_config=copy.deepcopy(self.variable_entries))
        #import pdb; pdb.set_trace()
        self.dual = self.create_dual()

    def create_dual(self):
        nodes_coords = self.primal.nodes.center[:]
        edges_center = self.primal.edges.center[:]
        faces_center = self.primal.faces.center[:]
        volumes_center = self.primal.volumes.center[:]
        sizes = np.array([len(nodes_coords), len(edges_center), len(faces_center), len(volumes_center)])
        sizes = np.cumsum(sizes)
        sizes = np.concatenate((np.array([0]), sizes))
        tag = lambda ind, type: (sizes[type] + ind +1)
        all_tetra = np.array([[], [], [], []]).T
        for vol in self.primal.volumes.all:
            faces = self.primal.volumes.adjacencies[vol]
            edges = self.primal.volumes._adjacencies(vol, dim_tag=1)
            nodes = self.primal.volumes.connectivities[vol]
            for edge in edges.T:
                adj_faces = np.intersect1d(self.primal.edges.bridge_adjacencies(edge, interface="edges",target="faces"), faces)
                node_edge = self.primal.edges.connectivities[edge.T].ravel()
                tetras = np.array([[tag(edge, 1), tag(node_edge[0], 0), tag(adj_faces[0], 2), tag(vol, 3)],
                                   [tag(edge, 1), tag(node_edge[0], 0), tag(adj_faces[1], 2), tag(vol, 3)],
                                   [tag(edge, 1), tag(node_edge[1], 0), tag(adj_faces[0], 2), tag(vol, 3)],
                                   [tag(edge, 1), tag(node_edge[1], 0), tag(adj_faces[1], 2), tag(vol, 3)]])
                all_tetra = np.vstack((all_tetra, tetras))
        print('Creating Coarse Scale Forming Dual Grid')
        dual = FineScaleMesh(mesh_file=None, dim=3, var_config=self.variable_entries)
        dual.core.mb.create_vertices(np.vstack((nodes_coords, edges_center, faces_center, volumes_center)))
        for tetra in all_tetra:
            dual.core.mb.create_element(types.MBTET, tetra.ravel().astype("uint64"))
        dual.core.run()
        dual.run()
        return dual


class simplePartition(object):
    def __init__(self, M):
        # num_of_vol, rx,ry,rz ,nx = 3, ny = 3, nz =3
        self.M = M

    def __call__(self, nx=4, ny=4, nz=4):
        return self.scheme(nx, ny, nz)

    def scheme(self, nx=4, ny=4, nz=4):
        centerCoord = self.M.volumes.center[:]
        num_of_vol = len(self.M)
        rx, ry, rz = self.M.rx, self.M.ry, self.M.rz
        if (rz[1] == 0) & (rz[0] == 0):
            nz = 1
            rz = (-1,1)
        box = np.array([0, (rx[1] - rx[0])/nx, 0,
              (ry[1] - ry[0]) /ny, 0,(rz[1]- rz[0])/(nz+0)]).reshape(3,2)
        cent_coord_El1 = box.sum(axis =1)/2
        tag = np.zeros(num_of_vol).astype("int")
        coarseCenters = np.zeros((nx*ny*nz,3))
        index = 0
        init_coords = np.array([rx[0],ry[0],rz[0]])
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    inc = np.multiply(box[:,1], np.array([x,y,z]))

                    #cent = cent_coord_El1 + inc
                    coarseCenters[index] = cent_coord_El1 + inc
                    # pdb.set_trace()

                    #inc = np.array([(nx) * x, (ny) * y, (nz) * z])
                    boxMin = box[:,0] + inc + init_coords
                    boxMax = box[:,1] + inc + init_coords
                    point = self.check_in_box(centerCoord,x=(boxMin[0], boxMax[0]), y=(boxMin[1], boxMax[1]) , z=(boxMin[2], boxMax[2]))
                    tag[point] = index
                    index += 1
        return self.tag_adjust(tag, coarseCenters)

    def check_in_box(self,coords, x , y, z):
        tag1 = (coords[:,0] > x[0]) & (coords[:,0] < x[1])
        tag2 = (coords[:,1] > y[0]) & (coords[:,1] < y[1])
        tag3 = (coords[:,2] > z[0]) & (coords[:,2] < z[1])
        return tag1 & tag2 & tag3

    def tag_adjust(self, tag, coarseCenter):
        fineTag =  tag
        elementsOriginal = [*set(tag)]
        elementsNovo = [*set(range(len(elementsOriginal)))]
        elementsMissing = set(range(len(coarseCenter))) - set(elementsOriginal)
        for elo, eln in zip(elementsOriginal,elementsNovo):
            if elo != eln:
                pointer = (tag == elo)
                fineTag[pointer] = eln
        return fineTag.astype(int) , np.delete(coarseCenter, [*elementsMissing], axis = 0)
