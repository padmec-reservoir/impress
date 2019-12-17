import numpy as np
cimport numpy as np

def hello():
    print('hello')

def point_in_volumes(vol_p_c, face_idx, p_coords, float tol = 0):
    cdef int i=0
    cdef int j=0
    cdef int k=0
    cdef int pts_count=0
    cdef bint outBox = False
    cdef int countOut=0
    cdef np.ndarray [np.float64_t, ndim = 2] vol_points_coords = vol_p_c
    cdef np.ndarray [np.int64_t, ndim = 2] faces_index = face_idx
    cdef np.ndarray [np.float64_t, ndim = 2] points_coords = p_coords
    cdef np.ndarray [np.float64_t, ndim = 2] lim = np.empty((3,2), dtype = np.float64) #linha: x/y/z #coluna: min/max
    cdef np.ndarray [np.uint8_t, cast = True, ndim = 1] retStatus = np.ones(dtype = np.uint8, shape = points_coords.shape[0])
    cdef np.ndarray [np.float64_t, ndim = 2] faces_centers = np.empty((faces_index.shape[0], 3), dtype = np.float64)
    cdef np.ndarray [np.float64_t, ndim = 1] average_point
    cdef np.ndarray [np.float64_t, ndim = 1] vol_center
    cdef np.ndarray [np.float64_t, ndim = 2] faces_normal = np.empty((faces_index.shape[0], 3),dtype = np.float64)
    cdef np.ndarray [np.float64_t, ndim = 1] v1_aux = np.empty(3, dtype=np.float64)
    cdef np.ndarray [np.float64_t, ndim = 1] v2_aux = np.empty(3, dtype=np.float64)
    #inicializa valores min/max
    for i in range(3):
      lim[i][0] = vol_points_coords[0][i]
      lim[i][1] = vol_points_coords[0][i]
    #encontra valores min/max
    for i in range(1,vol_points_coords.shape[0]):
      for j in range(3):
        if(vol_points_coords[i][j] < lim[j][0]):
          lim[j][0] = vol_points_coords[i][j]
        elif(vol_points_coords[i][j] > lim[j][1]):
          lim[j][1] = vol_points_coords[i][j]
    #checa pontos fora da boundingbox
    for i in range(points_coords.shape[0]):
      for j in range(3):
        outBox = outBox or (points_coords[i][j] < lim[j][0] or points_coords[i][j] > lim[j][1])
      if outBox:
        retStatus[i]=0
        countOut = countOut+1
        outBox=False
    #retorna todos os pontos estao fora
    if countOut == retStatus.size:
      print('All points are out of the boundingbox')
      return retStatus
    #calcula ponto medio do volume
    average_point = np.zeros(3, dtype = np.float64)
    for j in range(vol_points_coords.shape[0]):
      average_point[0] = average_point[0] + vol_points_coords[j][0]
      average_point[1] = average_point[1] + vol_points_coords[j][1]
      average_point[2] = average_point[2] + vol_points_coords[j][2]
    vol_center = average_point/vol_points_coords.shape[0]
    #calcula pontos medios e vetores normais das faces
    for i in range(faces_index.shape[0]):
      average_point = np.zeros(3, dtype = np.float64)
      #calculo do ponto medio
      pts_count = 0
      for j in range(faces_index[i].size):
        if faces_index[i][j]<0:
          break
        pts_count=pts_count+1
        average_point[0] = average_point[0] + vol_points_coords[faces_index[i][j]][0]
        average_point[1] = average_point[1] + vol_points_coords[faces_index[i][j]][1]
        average_point[2] = average_point[2] + vol_points_coords[faces_index[i][j]][2]
      average_point = average_point/pts_count
      for j in range(3):
        faces_centers[i][j] = average_point[j]
      #calculo do vetor normal
      for j in range(3):
        v1_aux[j]=vol_points_coords[faces_index[i][1]][j]-vol_points_coords[faces_index[i][0]][j]
        v2_aux[j]=vol_points_coords[faces_index[i][2]][j]-vol_points_coords[faces_index[i][0]][j]
      faces_normal[i]=np.cross(v1_aux, v2_aux)
      #checa se vetor normal esta apontando para fora
      for j in range(3):
        v1_aux[j]=faces_centers[i][j]-vol_center[j]
      if (v1_aux[0]*faces_normal[i][0]+v1_aux[1]*faces_normal[i][1]+v1_aux[2]*faces_normal[i][2] < 0 ):
        for j in range(3):
          faces_normal[i][j] = -faces_normal[i][j]

    #decide se cada ponto esta dentro (1), na face (2), ou fora (0)
    for i in range(points_coords.shape[0]):
      if(retStatus[i]):
        for j in range(faces_centers.shape[0]):
          for k in range(3):
            v1_aux[k] = points_coords[i][k]-faces_centers[j][k]
          prod = v1_aux[0]*faces_normal[j][0]+v1_aux[1]*faces_normal[j][1]+v1_aux[2]*faces_normal[j][2]
          if(prod > tol):
            retStatus[i]=0
            break
          elif(prod>=0):
            retStatus[i]=2
            break
    return retStatus
