

import numpy as np
import yaml
import os
from ..meshHandle.finescaleMesh import FineScaleMesh
# from ..meshHandle import imprutil as ip

# from ..CythonUtil import imprutil as ipq
from ..meshHandle.configTools.configClass import variableInit
from pymoab import types
import copy


def global_to_local(global_matrix):
    global_to_local_dic = {}
    global_connectivities_id = np.unique(global_matrix)
    local_connectivities_id = np.arange(len(global_connectivities_id))
    for key, val in zip(global_connectivities_id,
                        local_connectivities_id):
        global_to_local_dic[key] = val
    return np.vectorize(global_to_local_dic.get)(global_matrix)


def check_nodes_in_volume(el_coords, faces_connectivities, faces_normal, centers):
    el_center = el_coords.mean(axis=0).reshape((-1, 3))
    faces_center = el_coords[faces_connectivities].mean(axis=1)
    pseudo_normal = faces_center - el_center.repeat(len(faces_center), axis=0)
    change_sign = (faces_normal*pseudo_normal).sum(axis=1) <= 0
    faces_normal[change_sign] = -1 * faces_normal[change_sign]
    # look this up
    nodes_indicator = np.zeros((len(centers),
                                faces_connectivities.shape[0]), dtype=bool)
    for el in range(faces_connectivities.shape[0]):
        plane_check_flag = semi_plan_check(np.vstack((centers,el_center)), faces_normal[el],
                                           faces_center[el])
        if plane_check_flag[-1]:
            nodes_indicator[:, el] = plane_check_flag[:-1]
        else:
            nodes_indicator[:, el] = ~plane_check_flag[:-1]
    return nodes_indicator.all(axis=1)


def semi_plan_check(coords_list, normal_plane, point_on_plane, tol=1e-8):
    """coord_list -> list of points to be checked
    normal_plane -> normal vector to plane used for the comparision
    point_on_plane -> a point on the plane
    returns a bool of the size of coords_list with True or false """
    center_to_coords = coords_list - \
        np.repeat(point_on_plane.reshape((-1,3)), len(coords_list), axis=0)
    normal_plane = \
        np.repeat(normal_plane.reshape((-1,3)), len(coords_list), axis=0)
    inner_product = np.sum(center_to_coords*normal_plane,axis=1)
    flag = np.zeros(inner_product.shape, dtype=bool)

    #import pdb; pdb.set_trace()
    #flag[np.abs(inner_product) >= tol] = True
    flag[inner_product >= 0] = True
    # import pdb; pdb.set_trace()
    # flag[np.abs(inner_product) < tol] = True
    #flag[np.abs(inner_product) < tol] =True
    return flag


def check_in_box(coords, x, y, z):
    tag1 = (coords[:,0] >= x[0]) & (coords[:,0] <= x[1])
    tag2 = (coords[:,1] >= y[0]) & (coords[:,1] <= y[1])
    tag3 = (coords[:,2] >= z[0]) & (coords[:,2] <= z[1])
    return tag1 & tag2 & tag3


def tag_adjust(tag, coarseCenter):
    fineTag = tag
    elementsOriginal = [*set(tag)]
    elementsNovo = [*set(range(len(elementsOriginal)))]
    elementsMissing = set(range(len(coarseCenter))) - set(elementsOriginal)
    for elo, eln in zip(elementsOriginal,elementsNovo):
        if elo != eln:
            pointer = (tag == elo)
            fineTag[pointer] = eln
    return fineTag.astype(int), np.delete(coarseCenter, [*elementsMissing], axis = 0)


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
            self.arg = config_object.tree['Simple']['nx'], + \
                       config_object.tree['Simple']['ny'], + \
                       config_object.tree['Simple']['nz']
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
        self.variable_entries.add_entry('part', 'volumes', 1, 'int')
        self.print_creating('Primal Forming Mesh')
        self.primal = FineScaleMesh(file_path, var_config=copy.deepcopy(self.variable_entries))
        self.print_creating('Dual Forming Mesh')
        self.dual = self.create_forming_dual()
        volumes_indicator = self.find_primal_coarse_volumes()
        #import pdb; pdb.set_trace()
        #self.volumes_indicator_improvemnent(volumes_indicator)
        #import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()
        volumes_indicator = self.volumes_indicator_face_improv(volumes_indicator)
        # volumes_indicator = self.volumes_indicator_face_improv(volumes_indicator)
        # volumes_indicator = self.volumes_indicator_face_improv(volumes_indicator)
        #volumes_indicator = self.volumes_indicator_improvemnent(volumes_indicator)
        partition_tag = self.volumes_indicator_to_partition(volumes_indicator)
        proto_center = self.primal.volumes.center[:]
        # partition_tag = self.node_indicator_to_partition(node_indicator)
        # for el in range(node_indicator.shape[1]):
        #     self.M.nodein[node_indicator[:,el]] = 1*el
        # [proto_partition, proto_center] = self.find_primal_coarse_volumes()
        #import pdb; pdb.set_trace()
        #return tag_adjust(partition_tag, proto_center)
        return tag_adjust(partition_tag.ravel(), proto_center)

    def print_creating(self, text):
        print('-----------------------------------')
        print('  Creating {}'.format(text))
        print('-----------------------------------')
        pass

    def volumes_indicator_face_improv(self, volumes_indicator):
        all_adj = self.M.volumes.adjacencies[:]
        shared_volumes = volumes_indicator.sum(axis=1) > 1
        for coarse_el in range(volumes_indicator.shape[1]):
            analysed_elements = shared_volumes & volumes_indicator[:, coarse_el]
            local_adj, count = np.unique(all_adj[analysed_elements],
                                         return_counts=True)
            ref = 2*np.ones(len(self.M.faces))
            ref[local_adj] = count
            local_face_count = ref[all_adj[analysed_elements]]
            # import pdb; pdb.set_trace()
            # cond_1 = local_face_count == 1
            # cond_2 = local_face_count == 0
            # not_connected = (cond_1 | cond_2).all(axis=1)
            not_connected = (local_face_count.sum(axis=1) ==
                             local_face_count.shape[1])
            volumes_indicator[analysed_elements, coarse_el][not_connected]  = False
        return volumes_indicator

    def volumes_indicator_improvemnent(self, volumes_indicator):
        print(volumes_indicator)
        shared_volumes = volumes_indicator.sum(axis=1) > 1
        for index, el in enumerate(shared_volumes):
            if el:
                disputed_el_nodes = self.M.volumes.connectivities[index]
                disputing_coarse = volumes_indicator[index, :]
                disputing_coarse_tag = np.where(disputing_coarse)[0]
                volumes_in_tag = volumes_indicator[:, disputing_coarse]
                all_nodes = []
                for row in volumes_in_tag.T:
                    row_mod = np.copy(row)
                    row_mod[index] = False
                    shared_nodes = np.isin(disputed_el_nodes,np.unique(self.M.volumes.connectivities[row_mod.T])).sum()
                    all_nodes.append(shared_nodes)
                all_nodes = np.array(all_nodes)
                winner = disputing_coarse_tag[np.where(all_nodes == all_nodes.max())[0][0]]
                volumes_indicator[index, :] = False
                volumes_indicator[index, winner] = True
            print(index)
        return volumes_indicator
        # all_nodes.append(np.isin(disputed_el_nodes,np.unique(self.M.volumes.connectivities[row_mod.T]).sum()))

        # for coarse_el in range((volumes_indicator.shape[1])):
        #     # proto_partition = self.volumes_indicator_to_partition(volumes_indicator)
        #     check_volumes = volumes_indicator[:, coarse_el]
        #     el_multiple_coarse = volumes_indicator[check_volumes, :].sum(axis=1)
            # check_volumes[check_volumes] = (el_multiple_coarse > 1)
            # for fel in np.where(check_volumes)[0]:
            #     coarse_el_disputing = np.where(volumes_indicator[fel, :])[0]
            #     all_nodes = []
            #     for cel in coarse_el_disputing:
            #         all_nodes.append(np.unique(self.M.volumes.connectivities[volumes_indicator[:, cel]]))
            #     element_node = np.unique(self.M.volumes.connectivities[fel])
            #     shared_nodes = []
            #     for index in range(len(all_nodes)):
            #         import pdb; pdb.set_trace()
            #         shared_nodes.append(sum(np.isin(element_node, all_nodes[index])))
            #     import pdb; pdb.set_trace()

    def find_primal_coarse_volumes(self):
        all_fine_coords = self.M.nodes.coords[:]
        all_fine_centers = self.M.volumes.center[:]
        all_coords = self.primal.nodes.coords[:]
        connectivities = self.primal.volumes.connectivities[:]
        elements_coords = all_coords[connectivities]
        max_el = elements_coords.max(axis=1)
        min_el = elements_coords.min(axis=1)
        bbox_indicator = np.zeros((len(all_fine_coords),
                                   len(self.primal.volumes)), dtype=bool)
        cbox_indicator = np.zeros((len(all_fine_centers),
                                   len(self.primal.volumes)), dtype=bool)
        volumes_indicator = np.zeros((len(all_fine_centers),
                                      len(self.primal.volumes)), dtype=bool)
        x = (min_el[:, 0], max_el[:, 0])
        y = (min_el[:, 1], max_el[:, 1])
        z = (min_el[:, 2], max_el[:, 2])
        for el in range(len(self.primal.volumes)):
            x_el = (x[0][el], x[1][el])
            y_el = (y[0][el], y[1][el])
            z_el = (z[0][el], z[1][el])
            bbox_indicator[:, el] = check_in_box(all_fine_coords, x_el, y_el, z_el)
            cbox_indicator[:, el] = check_in_box(all_fine_centers, x_el, y_el, z_el)

        for el in range(len(max_el)):
            #import pdb; pdb.set_trace()
            el_node_ord = np.sort(self.primal.volumes.connectivities[el])
            local_coords_ord = self.primal.nodes.coords(el_node_ord)
            local_coords = elements_coords[el]
            faces_connectivities = \
                self.primal.faces.connectivities[self.primal.volumes.adjacencies[el]]
            faces_normal = self.primal.faces.normal[np.sort(self.primal.volumes.adjacencies[el])]
            local_faces_connectivities = global_to_local(faces_connectivities)
            cond_node = check_nodes_in_volume(local_coords_ord,
                                              local_faces_connectivities,
                                              faces_normal, all_fine_coords[bbox_indicator[:,el]])
            bbox_indicator[bbox_indicator[:, el], el] = cond_node
            cond_1 = bbox_indicator[:, el][self.M.volumes.connectivities[:]].any(axis=1)
            cond_2 = check_nodes_in_volume(local_coords_ord,
                                           local_faces_connectivities,
                                           faces_normal, all_fine_centers[cbox_indicator[:,el]])
            cbox_indicator[cbox_indicator[:, el], el] = cond_2
            cond_3 = np.logical_or(cond_1, cbox_indicator[:, el])
            volumes_indicator[:, el] = cond_3
        return volumes_indicator

    def volumes_indicator_to_partition(self, volumes_indicator):
        #import pdb; pdb.set_trace()
        partition_flag = np.zeros((volumes_indicator.shape[0], 1), dtype=int)
        #
        # for el in range(len(self.primal.volumes)):
        #     partition_indicator[:,el] = node_indicator[fine_connectivities,el].any(axis=1)
        for el in range(len(self.primal.volumes)):
            partition_flag[volumes_indicator[:, el]] = 1*el
        #import pdb; pdb.set_trace()
        return partition_flag

    # def find_primal_coarse_volumes3(self):
    #     all_coords = self.primal.nodes.coords[:]
    #     connectivities = self.primal.volumes.connectivities[:]
    #     elements_coords = all_coords[connectivities]
    #     max_el = elements_coords.max(axis=1)
    #     min_el = elements_coords.min(axis=1)
    #     fine_scale_elements_coords = self.M.volumes.center[:]
    #     bbox_indicator = np.zeros((len(fine_scale_elements_coords),
    #                                len(self.primal.volumes)), dtype=bool)
    #     coarse_indicator = np.zeros((len(fine_scale_elements_coords),
    #                                len(self.primal.volumes)), dtype=bool)
    #     primal_indicator = np.zeros((len(fine_scale_elements_coords), 1), dtype=int)
    #     for el in range(len(max_el)):
    #         x = (min_el[el, 0], max_el[el, 0])
    #         y = (min_el[el, 1], max_el[el, 1])
    #         z = (min_el[el, 2], max_el[el, 2])
    #         bbox_indicator[:, el] = check_in_box(fine_scale_elements_coords,
    #                                              x, y, z)
    #     index = 0
    #     for el in range(len(max_el)):
    #         local_coords = elements_coords[el]
    #         faces_connectivities = \
    #             self.primal.faces.connectivities[self.primal.volumes.adjacencies[el]]
    #
    #         # faces_connectivities = \
    #         #     1+1
    #         # faces_edges = self.primal.faces.adjacencies[self.primal.volumes.adjacencies[el]]
    #         # local_faces_edges = global_to_local(faces_edges)
    #         # local_inner_edges = global_to_local(self.primal.edges.adjacencies[np.unique(faces_edges)])
    #         # local_inner_edges[local_faces_edges]
    #         faces_normals = self.primal.faces.normal[self.primal.volumes.adjacencies[el]]
    #         local_faces_connectivities = global_to_local(faces_connectivities)
    #         # local_faces_connectivities = \
    #         #     np.vectorize(global_to_local_dic.get)(faces_connectivities)
    #         # #import pdb; pdb.set_trace()
    #         in_volumes = ipq.point_in_volumes(local_coords,
    #                                           local_faces_connectivities,
    #                                           self.M.volumes.center[bbox_indicator[:,el]],0)
    #         # in_volumes = ipq.point_in_volumes(local_coords,
    #         #                                   local_faces_connectivities,
    #         #                                   self.M.volumes.center[:],0)
    #         #import pdb; pdb.set_trace()
    #         #import pdb; pdb.set_trace()
    #         #in_volumes2 = check_in_volume(local_coords, local_faces_connectivities, faces_normals, self.M.volumes.center[:])
    #         coarse_indicator[bbox_indicator[:, el], el] = (in_volumes >= 1)
    #         #import pdb; pdb.set_trace()
    #         # coarse_indicator[:,el] =  (in_volumes >= 1)
    #         #indices = np.where(bbox_indicator[:, el])[0][in_volumes != 0]
    #         primal_indicator[bbox_indicator[:, el]] = (el * (in_volumes >= 1)).reshape((-1,1))
    #         #import pdb; pdb.set_trace()
    #         #import pdb; pdb.set_trace()
    #         index += 1
    #     return primal_indicator, self.primal.volumes.center[:]

    def create_forming_dual(self):
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
                    coarseCenters[index] = cent_coord_El1 + inc
                    boxMin = box[:,0] + inc + init_coords
                    boxMax = box[:,1] + inc + init_coords
                    point = check_in_box(centerCoord,x=(boxMin[0], boxMax[0]),
                                         y=(boxMin[1], boxMax[1]), z=(boxMin[2]
                                         , boxMax[2]))
                    tag[point] = index
                    index += 1
        return tag_adjust(tag, coarseCenters)
