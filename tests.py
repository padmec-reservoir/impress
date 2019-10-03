import numpy as np
import matplotlib.pyplot as plt

import networkx as nx

x = 15

G = nx.graph.Graph()


neighbors = np.setdiff1d(M.coarse.elements[x].nodes.coarse_neighbors, len(M.coarse.elements))
totaledges = np.array([])
for el in neighbors:
    tag1 = M.coarse.iface_coarse_neighbors[:,0] == el
    tag2 = np.isin(M.coarse.iface_coarse_neighbors[:,1],neighbors)
    tag3 = M.coarse.iface_coarse_neighbors[:,0] != M.coarse.iface_coarse_neighbors[:,1]
    tag4 = (M.coarse.iface_coarse_neighbors[:,0] != x) | (M.coarse.iface_coarse_neighbors[:,1] != x)
    indices = np.where(tag1 & tag2 & tag3 & tag4)[0]
    #print(M.coarse.iface_coarse_neighbors[indices,:])
    L = M.coarse.iface_coarse_neighbors[indices, :]
    #import pdb; pdb.set_trace()
    for line in L:
        print(line)
        #import pdb; pdb.set_trace()
        # if line[0] is not line[1]:
        G.add_edge(line[0], line[1])
    for index in indices:
        #import pdb; pdb.set_trace()
        totaledges = np.append(totaledges, l.coarse_edges[index].ravel())


M.dreams[:] = -1
M.dreams[np.unique(totaledges)] = 2
M.dreams[l.coarse_center[neighbors]] = 3
M.core.print(folder = "results")
nx.draw(G,with_labels=True)
plt.show()
import networkx.algorithms as alg
    #import pdb; pdb.set_trace()
    # part = M.coarse.partition[:].ravel()[l.interface_vol]
# for index in range(M.coarse.num_internal_faces):
#     line = part[index]
#     G.add_edge(line[0], line[1])

# import numpy as np
# import pdb
# M.dreams[:] = -1
# # for index, el in enumerate(l.coarse_edges):
# #     print(index,el)
# #     M.dreams[el] = index
#
# # for index, el in enumerate(l.coarse_edges):
# #     M.dreams[el] = index
# num = 80
#
# for index, el in enumerate(l.coarse_edges):
#     #pdb.set_trace()
#     if index == num:
#         M.dreams[el] = np.arange(len(el)).astype(float).T
#
# M.core.print()
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
