import numpy as np
import matplotlib.pyplot as plt

import networkx as nx
#x = 10
x = 13

import networkx.algorithms as alg
G = nx.graph.Graph()
neighbors = np.setdiff1d(M.coarse.elements[x].nodes.coarse_neighbors, len(M.coarse.elements))
totaledges = np.array([])
print("Neighbors of Coarse Volume: " + str(x))
print(neighbors)
for el in neighbors:
    print("Planar Space of Neighbors of Coarse Volume " + str(el))
    planar_local_neighbors = np.intersect1d(M.coarse.elements[el].nodes.coarse_neighbors,neighbors)
    print(planar_local_neighbors)
    G = nx.graph.Graph()

    for element in planar_local_neighbors:

        tag1 = M.coarse.iface_coarse_neighbors[:,0] == element
        tag2 = np.isin(M.coarse.iface_coarse_neighbors[:,1],planar_local_neighbors)
        tag3 = M.coarse.iface_coarse_neighbors[:,0] != M.coarse.iface_coarse_neighbors[:,1]
        tag4 = (M.coarse.iface_coarse_neighbors[:,0] != el) | (M.coarse.iface_coarse_neighbors[:,1] != el)
        indices = np.where(tag1 & tag2 & tag3 & tag4)[0]
        L = M.coarse.iface_coarse_neighbors[indices, :]
    #import pdb; pdb.set_trace()
        for line in L:
            #print(line)
            #import pdb; pdb.set_trace()
            # if line[0] is not line[1]:
            G.add_edge(line[0], line[1])
        for index in indices:
            #import pdb; pdb.set_trace()
            totaledges = np.append(totaledges, l.coarse_edges[index].ravel())
    print(alg.planarity.check_planarity(G))
    # if not alg.planarity.check_planarity(G)[0]:
    #     import pdb; pdb.set_trace()

    nx.draw(G,with_labels=True)
    plt.show()
# #
# M.dreams[:] = -1
# M.dreams[np.unique(totaledges)] = 2
# M.dreams[l.coarse_center[neighbors]] = 3
M.core.print(folder = "results")
# nx.draw(G,with_labels=True)
# plt.show()
# import networkx.algorithms as alg
