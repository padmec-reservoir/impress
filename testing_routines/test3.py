from scipy.spatial.qhull import Delaunay

import numpy as np
import matplotlib.pyplot as plt
import tools.view as tl
import networkx as nx

x = 16

G = nx.graph.Graph()


#finding neighbors of the x element
neighbors = np.setdiff1d(M.coarse.elements[x].nodes.coarse_neighbors, len(M.coarse.elements))
neighbors = np.sort(np.append(neighbors,x)).astype(int)

#variable containing which elements are on the boundary
boundary_test_neighbors = np.zeros(neighbors.shape, dtype=bool)
for index, el in enumerate(neighbors):
    boundary_test_neighbors [index] = M.coarse.elements[index].faces.is_on_father_boundary


faces_neighbors = np.isin(M.coarse.partition[:].ravel()[l.interface_vol], neighbors)
#import pdb; pdb.set_trace()

boundary_ifaces = np.arange(len(M.coarse.interfaces_faces)) > M.coarse.num_internal_faces

faces_neighbors_internal = np.where(faces_neighbors[:, 0] & faces_neighbors[:, 1] & ~ boundary_ifaces)[0]
#import pdb; pdb.set_trace()
faces_neighbors_external = np.where(faces_neighbors[:, 0] & faces_neighbors[:, 1] & boundary_ifaces)[0]
extra_coords_internal = M.faces.center[l.center_face[faces_neighbors_external]]
extra_coords_external = M.faces.center[l.center_face[faces_neighbors_internal]]

#import pdb; pdb.set_trace()
neighbors_ref = np.arange(0,len(neighbors))
neighbors_ref = np.append(neighbors_ref, -2)
neighbors_ref = np.append(neighbors_ref, -1)
#import pdb; pdb.set_trace()
neighbors_ext = np.hstack((neighbors, -1 * np.ones(len(extra_coords_internal)), -2 * np.ones(len(extra_coords_external))))


#creating triangulation coordinates
dcoords = M.volumes.center[l.coarse_center[neighbors]]
ncoords = np.vstack((dcoords, extra_coords_internal, extra_coords_external))

#perfoming delaunay triangulation
trian = Delaunay(ncoords)

#find boundary nodes
boundary_nodes = np.unique(trian.convex_hull)

#finding boundary tetrahedron - tetrahedron that share one face with the boundary_edges
boundary = np.isin(trian.simplices, boundary_nodes)
tet_on_boundary = np.sum(boundary,axis = 1) >= 3
btet = neighbors_ext[trian.simplices[tet_on_boundary]]
tag1 = np.sum(btet == -2, axis=1) >=1
tag2 = np.sum(btet == -1, axis=1) == 3
tag3 = np.sum(btet == -1, axis=1) == 4
remove1 = tag2 & tag1
remove2 = tag3
remove_index = np.unique(np.hstack((np.where(tet_on_boundary)[0][remove1], np.where(tet_on_boundary)[0][remove2])))


import pdb; pdb.set_trace()

# missing = np.setdiff1d(np.setdiff1d(neighbors_ref, boundary_nodes), np.where(neighbors == x)[0])
# # finding tetrahedrons that share at least a face with the boundary
# boundary_tetra = np.sum(np.isin(trian.simplices, boundary_nodes), axis=1) >= 3

#bad_tetra = np.zeros(boundary_tetra.shape, dtype=bool)

index = 0

# for index, value in enumerate(boundary_tetra):
#     if value:
#         #print('Corno')
#         local_el = trian.simplices[index]
#         # import pdb; pdb.set_trace()
#         local_missing = np.setdiff1d(np.setdiff1d(local_el, boundary_nodes),np.where(neighbors == x)[0])
#         if np.isin(local_missing, missing):
#             bad_tetra[index] = True
#             #break

# calculating volume of the tetrahedrons
volumes = np.zeros(len(trian.simplices))
index = 0
for line in trian.simplices:
    point = trian.points[line]
    volumes[index] = np.abs(tl.volume(point[0], point[1], point[2], point[3]))
    index += 1
    #import pdb; pdb.set_trace()

trian.simplices = np.delete(trian.simplices, remove_index, axis=0)
#good_tetra = trian.simplices
        #print(index, value)
        #print(trian.simplices[index])
    #print(index,line)

# for line, ref in zip(trian.simplices, ~np.isin(trian.simplices[boundary_tetra], boundary_nodes)) :
#     print(line,ref, line[ref])
#     if np.isin(line[ref], missing):
#         bad_tetra[index] = True
#     index += 1
a = tl.DelaunayView(ncoords, good_tetra, neighbors, x)

# xyz = dcoords
# import mpl_toolkits.mplot3d as a3
# axes = a3.Axes3D(pl.figure())
# vts = xyz[tri, :
#
# nodes = 1

# totaledges = np.array([])
# for el in neighbors:
#     tag1 = M.coarse.iface_coarse_neighbors[:,0] == el
#     tag2 = np.isin(M.coarse.iface_coarse_neighbors[:,1],neighbors)
#     tag3 = M.coarse.iface_coarse_neighbors[:,0] != M.coarse.iface_coarse_neighbors[:,1]
#     tag4 = (M.coarse.iface_coarse_neighbors[:,0] != x) | (M.coarse.iface_coarse_neighbors[:,1] != x)
#     indices = np.where(tag1 & tag2 & tag3 & tag4)[0]
#     #print(M.coarse.iface_coarse_neighbors[indices,:])
#     L = M.coarse.iface_coarse_neighbors[indices, :]
#     #import pdb; pdb.set_trace()
#     for line in L:
#         print(line)
#         #import pdb; pdb.set_trace()
#         # if line[0] is not line[1]:
#         G.add_edge(line[0], line[1])
#     for index in indices:
#         #import pdb; pdb.set_trace()
#         totaledges = np.append(totaledges, l.coarse_edges[index].ravel())
#
#
# M.dreams[:] = -1
# M.dreams[np.unique(totaledges)] = 2
# M.dreams[l.coarse_center[neighbors]] = 3
# M.core.print(folder = "results")
# nx.draw(G,with_labels=True)
# plt.show()
# import networkx.algorithms as alg
