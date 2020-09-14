"""
Schemes of partitioning
"""

#import pdb
import numpy as np
from numba import jit

#para definir funções de particionamento deve-se atentar as saidas
#é necessário retornar um vetor tag que associa uma partição a cada celula fina
#também é necessario retornar o uma matrix
# n x 3 com o centro de cada volume coarse, n -> numero de elementos do volume coarse

#a função tagAdjust remove renumera os elementos da malha coarse de maneira que eles permanecam
# continuos contando de 0 ao N de volumes coarse, a função tb remover os centros indesejados

#para adicionar novos esquemas basta criar um esquema com numeracao sequencial
#ex. scheme2  , mais shceme3
#a seção campo de leitura do arquivo msCoarse.ini deverá ter nome correspondente
#ex [Coarsening_2_Input] e [Coarsening_3_Input]
#todos atributos serao passados para funcao correspondente na ordem definida
#por default a leitura dos elementos é float, caso necessário. converta para int
#ex int(nx)

def scheme1(centerCoord, num_of_vol, rx,ry,rz ,nx = 3, ny = 3, nz =3 ):
    #input : centerCoord - > array with the center of elements
    #        num_of_vol = number of volumes

    #        rx,ry,rz - (min, max) values of x,y,z of the phyisical domain
    #        nx, ny, nz
    # msh -> objeto da clase meshUtil
    #centerCoord = msh.readData("CENTER")
    nx = int(nx)
    ny = int(ny)
    nz = int(nz)
    if (rz[1] == 0)  &  (rz[0] == 0):
        nz = 1
        rz = (-1,1)
    box = np.array([0, (rx[1] - rx[0])/nx, 0,(ry[1] - ry[0]) /ny, 0,(rz[1] - rz[0])/(nz+0)]).reshape(3,2)
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
                point = checkinBox(centerCoord,x=(boxMin[0], boxMax[0]), y=(boxMin[1], boxMax[1]) , z=(boxMin[2], boxMax[2]))
                tag[point] = index
                index += 1
    return tagAdjust(tag,coarseCenters)


def scheme2(centerCoord, num_of_vol, rx,ry,rz ,nx = 3, ny = 3, nz =3 ):
    #input : centerCoord - > array with the center of elements
    #        num_of_vol = number of volumes

    #        rx,ry,rz - (min, max) values of x,y,z of the phyisical domain
    #        nx, ny, nz
    # msh -> objeto da clase meshUtil
    #centerCoord = msh.readData("CENTER")
    print("ESQUEMA TIPO 2 PARTICIONAMENTO")
    nx = int(nx)
    ny = int(ny)
    nz = int(nz)
    box = np.array([0, (rx[1] - rx[0])/nx, 0,(ry[1] - ry[0]) /ny, 0,(rz[1] - rz[0])/(nz+0)]).reshape(3,2)
    cent_coord_El1 = box.sum(axis =1)/2
    tag = np.zeros(num_of_vol).astype("int")
    coarseCenters = np.zeros((nx*ny*nz,3))
    index = 0

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                inc = np.multiply(box[:,1], np.array([x,y,z]))
                #cent = cent_coord_El1 + inc
                coarseCenters[index] = cent_coord_El1 + inc
                #inc = np.array([(nx) * x, (ny) * y, (nz) * z])
                boxMin = box[:,0] + inc
                boxMax = box[:,1] + inc
                point = checkinBox(centerCoord,x=(boxMin[0], boxMax[0]), y=(boxMin[1], boxMax[1]) , z=(boxMin[2], boxMax[2]))
                tag[point] = index
                index += 1
    #return tagAdjust(tag,coarseCenters)


def scheme3(centerCoord, num_of_vol, rx,ry,rz ,nx = 3, ny = 3, nz =3 ):
    #input : centerCoord - > array with the center of elements
    #        num_of_vol = number of volumes

    #        rx,ry,rz - (min, max) values of x,y,z of the phyisical domain
    #        nx, ny, nz
    # msh -> objeto da clase meshUtil
    #centerCoord = msh.readData("CENTER")
    print("ESQUEMA 3")
    nx = int(nx)
    ny = int(ny)
    nz = int(nz)
    box = np.array([0, (rx[1] - rx[0])/nx, 0,(ry[1] - ry[0]) /ny, 0,(rz[1] - rz[0])/(nz+0)]).reshape(3,2)
    cent_coord_El1 = box.sum(axis =1)/2
    tag = np.zeros(num_of_vol).astype("int")
    coarseCenters = np.zeros((nx*ny*nz,3))
    index = 0
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                inc = np.multiply(box[:,1], np.array([x,y,z]))
                #cent = cent_coord_El1 + inc
                coarseCenters[index] = cent_coord_El1 + inc
                #inc = np.array([(nx) * x, (ny) * y, (nz) * z])
                boxMin = box[:,0] + inc
                boxMax = box[:,1] + inc
                point = checkinBox(centerCoord,x=(boxMin[0], boxMax[0]), y=(boxMin[1], boxMax[1]) , z=(boxMin[2], boxMax[2]))
                tag[point] = index
                index += 1
    return tagAdjust(tag,coarseCenters)



#checa se um ponto esta dentro de um cubo

@jit(parallel=True)
def checkinBox(coords, x , y, z):
    tag1 = (coords[:,0] > x[0])   &  (coords[:,0] < x[1])
    tag2 = (coords[:,1] > y[0])   &  (coords[:,1] < y[1])
    tag3 = (coords[:,2] > z[0])   &  (coords[:,2] < z[1])
    return tag1 & tag2 & tag3



#função para corrigir os tags em caso que as malhas geradoras possuam volumes grossos sem celulas dentro
#e remover as respectivas coordenadas do centro dos volumes das malhas primais
#  @jit(parallel=True, cache=True)
def tagAdjust(tag, coarseCenter):
    # msh -> objeto da clase meshUtil
    fineTag =  tag
    elementsOriginal = [*set(tag)]
    elementsNovo = [*set(range(len(elementsOriginal)))]
    elementsMissing = set(range(len(coarseCenter))) - set(elementsOriginal)
    for elo, eln in zip(elementsOriginal,elementsNovo):
        if elo != eln:
            pointer = (tag == elo)
            fineTag[pointer] = eln
    return fineTag.astype(int) , np.delete(coarseCenter, [*elementsMissing], axis = 0)
