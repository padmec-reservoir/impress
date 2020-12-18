from scipy.spatial.qhull import Delaunay
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tools.view as tl
import networkx as nx
from scipy.spatial import ConvexHull
x = 22

G = nx.graph.Graph()
#finding neighbors of the x element
neighbors = np.setdiff1d(M.coarse.elements[x].nodes.coarse_neighbors, len(M.coarse.elements))
neighbors = np.sort(np.insert(neighbors,1,x)) # XXX: ))
#neighbors = np.append(neighbors, x).astype(int)
#center coordinates of the CV
ccoord = M.volumes.center[l.coarse_center[x]]
#center coordinates of the cv coarse volume neigbhors
dcoords = M.volumes.center[l.coarse_center[np.setdiff1d(neighbors,x)]]
faces_neighbors = np.isin(M.coarse.partition[:].ravel()[l.interface_vol], neighbors)
all_neighbors_index = np.logical_and(faces_neighbors[:,0], faces_neighbors[:,1])
internal_tag = np.zeros(all_neighbors_index.shape, dtype="bool")
internal_tag[0:M.coarse.num_internal_faces] = 1

# removing bad tetrahedron that do not have the center
#modified_simplices = trian.simplices[np.sum(trian.simplices == trian.simplices.max(),axis=1) == 1]



#coordinates of the interfaces of the CV with its neighbors
x_tag = np.sum((M.coarse.partition[:].ravel()[l.interface_vol] == x),axis=1) >=1
cv_interfaces = np.logical_and(all_neighbors_index,x_tag)
cv_icoords = M.faces.center[l.center_face[np.logical_and(cv_interfaces, internal_tag)]]
#coordinates of the interface of the CV with the boundary
cv_xcoords = M.faces.center[l.center_face[np.logical_and(cv_interfaces, ~internal_tag)]]
#cv_neighbors = np.sum(np.isin(M.coarse.partition[:].ravel()[l.interface_vol],x),axis=1) >= 1


#coordinates of the interfaces of the CV neighbors among each other
icoords = M.faces.center[l.center_face[np.logical_and(~cv_interfaces, internal_tag)]]
#coordinates of the interfaces of the CV neighbors on the boundary
xcoords = M.faces.center[l.center_face[np.logical_and(~cv_interfaces, ~internal_tag)]]

#
# cv_interfaces = np.logical_and(all_neighbors_index, cv_neighbors)
# neighbors_interfaces = np.logical_and(all_neighbors_index, ~cv_neighbors)
#
# np.sum(np.isin(M.coarse.partition[:].ravel()[l.interface_vol],x),axis=1) == 1
#
# #coordinates of the interfaces of the CV with its neighbors
# cv_icoords = M.faces.center[l.center_face[np.logical_and(cv_interfaces, internal_tag)]]
# #coordinates of the interface of the CV with the boundary
# cv_xcoords = M.faces.center[l.center_face[np.logical_and(cv_interfaces, ~internal_tag)]]
#
# #coordinates of the interfaces of the CV neighbors among each other
# icoords = M.faces.center[l.center_face[np.logical_and(neighbors_interfaces, internal_tag)]]
# #coordinates of the interfaces of the CV neighbors on the boundary
# xcoords = M.faces.center[l.center_face[np.logical_and(neighbors_interfaces, ~internal_tag)]]

tcoords = np.vstack((dcoords, cv_icoords, cv_xcoords, icoords, xcoords, ccoord))
point_type = np.concatenate((0*np.ones(dcoords.shape[0]),1*np.ones(cv_icoords.shape[0]), 2*np.ones(cv_xcoords.shape[0]) , 3*np.ones(icoords.shape[0]), 4*np.ones(xcoords.shape[0]), 5*np.ones(ccoord.shape[0])))


tcoords = np.vstack((dcoords, cv_icoords, cv_xcoords, icoords, ccoord))
point_type = np.concatenate((0*np.ones(dcoords.shape[0]),1*np.ones(cv_icoords.shape[0]), 2*np.ones(cv_xcoords.shape[0]) , 3*np.ones(icoords.shape[0]),  5*np.ones(ccoord.shape[0])))


tcoords = np.vstack((dcoords, ccoord))
point_type = np.concatenate((0*np.ones(dcoords.shape[0]), 5*np.ones(ccoord.shape[0])))





trian = Delaunay(tcoords)

df = pd.DataFrame(point_type[trian.simplices])
df.to_excel("delaunay.xlsx", index=False)
modified_simplices = trian.simplices[np.sum(trian.simplices == trian.simplices.max(),axis=1) == 1]
a = tl.DelaunaySingle(trian.points, trian.simplices, point_type, 2)
#a = tl.DelaunaySingle(trian.points, trian.simplices,point_type,2)
# import pdb; pdb.set_trace()
#ncoords = np.vstack((dcoords, extra_coords_internal, extra_coords_external))

# a = tl.DelaunaySingle(, good_tetra, neighbors, x)
# # #perfoming delaunay triangulation
# # trian = Delaunay(dcoords)
