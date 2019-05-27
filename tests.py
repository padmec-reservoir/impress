import numpy as np
import pdb
M.dreams[:] = -1
# for index, el in enumerate(l.coarse_edges):
#     print(index,el)
#     M.dreams[el] = index

# for index, el in enumerate(l.coarse_edges):
#     M.dreams[el] = index
num = 23

for index, el in enumerate(l.coarse_edges):
    #pdb.set_trace()
    if index == num:
        M.dreams[el] = np.arange(len(el)).astype(float).T

M.core.print()
# from numba import jit
# import numpy as np
#
#
#
# def teste(qarray, vector):
#     print(np.asarray(qarray)[vector])
#
# n = M.core.all_volumes
# ll = np.array([1,3,4, 10, 35,37,99]); start = time.time()
# for i in range(1000):
#     teste(n, ll)
# end = time.time()
# print(start-end)
# M.happiness[:] = -1
# for x in range(len(M.coarse.interfaces_nodes)):
#     M.happiness[M.coarse.interfaces_nodes[x]] = x
#
#
#
# M.joy[:] = -1
# for x in range(len(M.coarse.interfaces_edges)):
#     M.joy[M.coarse.interfaces_edges[x]] = x
#
#
# M.pride[:] = -1
# for x in range(len(M.coarse.interfaces_faces)):
#     M.pride[M.coarse.interfaces_faces[x]] = x
#
#
# M.dreams[:] = -1
# for x in range(len(M.coarse.interfaces_faces)):
#     M.pride[M.coarse.interfaces_faces[x]] = x
