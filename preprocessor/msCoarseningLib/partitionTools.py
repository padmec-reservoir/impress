

import numpy as np
import yaml
import os
from ..meshHandle.finescaleMesh import FineScaleMesh


class configManager(object):
    def __init__(self, empty=False, file_path=None):
        if file_path == None:
            file_path = 'input_cards/coarsening.yml'
        if empty == False:
            self.read(file_path)
        else:
            self.tree = {'Scheme': '', 'Simple': '', 'Smart': ''}

    def scheme(self, scheme = 'smart'):
        self.tree['Scheme'] =  scheme

    def simple(self, nx=4, ny=4, nz=4):
        self.tree['Simple'] = {'nx': nx, 'ny': ny, 'nz': nz}

    def smart(self, path='default', file='mesh.h5m'):
        self.tree['Smart'] = {'path': path, 'file': file}

    def read(self, file_path):
        with open(file_path) as file: # Use file to refer to the file object
            data = file.read()
            self.tree = yaml.safe_load(data)

    def write(self, file_path):
        with open(file_path, 'w') as f:
                data = yaml.dump(self.tree, f)


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
        else:
            print('Scheme ' + scheme + ' not defined.')

    def __call__(self):
        return self.run()

    def run(self, *args):
        if len(args) == 0:
            return self.func(*self.arg)
        else:
            return self.func(*args)

    def smart(self, file_path):
        return self.partitioner(file_path)

    def simple(self, nx, ny, nz):
        return self.partitioner(nx, ny, nz)


class smartPartition(object):
    def __init__(self,M):
        self.M = M

    def __call__(self, file_path):
        return self.scheme(file_path)

    def scheme(self, file_path):
        self.primal = FineScaleMesh(file_path)
        self.pcenter = self.primal.volumes.center[:]
        
        import pdb; pdb.set_trace()




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
        tag1 = (coords[:,0] > x[0])   &  (coords[:,0] < x[1])
        tag2 = (coords[:,1] > y[0])   &  (coords[:,1] < y[1])
        tag3 = (coords[:,2] > z[0])   &  (coords[:,2] < z[1])
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
