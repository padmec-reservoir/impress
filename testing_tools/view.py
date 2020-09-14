from pymoab import core, types, rng, topo_util
from pymoab import skinner as sk
import numpy as np


def volume(p1, p2, p3, p4):
    d1 = p2 - p1
    d2 = p3 - p1
    d3 = p4 - p1
    matrix = np.vstack((d1, d2, d3))
    return np.abs((1/6)*np.linalg.det(matrix))


def volume_p(tet):
    vol = np.zeros(tet.shape[0])
    #import pdb; pdb.set_trace()
    for i in range(tet.shape[0]):
        p1, p2, p3, p4 = tet[i][0], tet[i][1], tet[i][2], tet[i][3]
        # d1 = p2 - p1
        # d2 = p3 - p1
        # d3 = p4 - p1
        matrix = np.hstack((np.vstack((p1, p2, p3, p4)), np.ones((4, 1))))
        vol[i] = (1/6)*np.linalg.det(matrix)
    return np.abs(vol)
    #matrix = np.vstack((d1, d2, d3))
    #return (1/6)*np.linalg.det(matrix)


class DelaunaySingle(object):
    def __init__(self, coords, elements, point_type, center):
        self.mb = core.Core()
        #self.point_type = np.concatenate((np.array([-1]), point_type))
        self.point_type = point_type
        verts = self.mb.create_vertices(coords)
        self.rs = self.mb.get_root_set()

        elements = elements + np.ones(elements.shape)
        self.tag_handle = self.mb.tag_get_handle("Teste", 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, create_if_missing = True)
        self.tag_node = self.mb.tag_get_handle("Nodes", 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, create_if_missing = True)

        for el in elements:
            tetra = self.mb.create_element(types.MBTET, el.ravel().astype("uint64"))

        # self.mtu = topo_util.MeshTopoUtil(self.mb)
        # self.mtu.construct_aentities(verts)
        #self.mtu.construct_aentities(verts)
        #self.elements = self.mb.get_entities_by_dimension(0, 3)
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        self.all_elements = self.mb.get_entities_by_dimension(0, 3)
        self.nodes = self.mb.get_entities_by_dimension(0, 0)
        self.faces = self.skinner(self.all_elements)
        self.boundary_nodes = self.find_boundary_nodes()
        self.boundary_elements = self.find_boundary_elements()

        self.desired_bnodes = np.array(self.nodes)[self.point_type == 0]
        self.find_simple_bad()
        #self.select(self.boundary_nodes)
        #self.remove_vol(self.all_elements[0:120])
        #self.find_bad_vol()
        #import pdb; pdb.set_trace()
        #self.smooth()
        #import pdb; pdb.set_trace()
        self.create_vtk()

    def find_missing_nodes(self):
        return self.desired_bnodes[~np.isin(self.desired_bnodes,np.array(self.boundary_nodes))]

    def bad_vol(self):
        index = 0
        point_type = np.concatenate((np.array([-1]), self.point_type))
        while not self.condition():
            #import pdb; pdb.set_trace()
            index = index + 1
            connect = self.connectivities(self.boundary_elements, 4)
            flag = point_type[connect]
            all_tetra_ref = np.arange(connect.shape[0])
            num_tetra_bnodes = np.sum(np.isin(connect,np.array(self.boundary_nodes)), axis=1)
            four_bnodes_tetra = (num_tetra_bnodes == 4)
            three_bnodes_tetra = (num_tetra_bnodes == 3)


            # treating tetrahedron with 4 nodes on the boundary those tetrahedron
            # will be removed as no node are exclusive to these entities
            rm1 = all_tetra_ref[four_bnodes_tetra][~np.any(np.logical_and(connect[four_bnodes_tetra],flag[four_bnodes_tetra] == 0),axis=1)]
            #~np.any(np.logical_and(connect[four_bnodes_tetra],flag[four_bnodes_tetra] == 0),axis=1)
            #import pdb; pdb.set_trace()

            # treating tetrahedron with 3 nodes on the boundary those tetrahedron
            # will be removed as no node are exclusive to these entities
            rmp = all_tetra_ref[three_bnodes_tetra]
            # ~np.isin(connect[three_bnodes_tetra],np.array(self.boundary_nodes))
            tag = np.any(flag[three_bnodes_tetra] == 0,  axis=1)
            rm2 = rmp[~tag]
            tetra_tag = rmp[tag]
            rm3 = tetra_tag[np.all(~np.isin(connect[tetra_tag],np.array(self.boundary_nodes)) == (flag[tetra_tag] == 0), axis=1)]
            # flag[three_bnodes_tetra]
            # rm2 = all_tetra_ref[four_bnodes_tetra]
            remove_tetra = np.sort(np.hstack((rm1,rm2,rm3)))
            remove_flag = np.isin(all_tetra_ref, remove_tetra)

            self.remove_vol(self.select(self.boundary_elements,remove_flag))

            self.create_vtk()
            #import pdb; pdb.set_trace()

    def find_simple_bad(self):
        index = 0
        point_type = np.concatenate((np.array([-1]), self.point_type))
        missing = self.find_missing_nodes()
        coord = self.mb.get_coords(self.nodes).reshape(-1,3)]
        remov = rng.Range()
        while not self.condition():
            index = index + 1
            connect = self.connectivities(self.boundary_elements, 4)
            all_tetra_ref = np.arange(connect.shape[0])
            flag = point_type[connect]
            for el in missing:
                tetra_analize = np.any(np.isin(connect, el), axis=1)
                np.unique(connect[tetra_analize])
                tet_coord = (connect[tetra_analize] - np.ones(connect[tetra_analize].shape)).astype("uint64")
                #tet_coord1 = (connect - np.ones(connect.shape)).astype("uint64")
                vol_tetra = volume_p(coord[tet_coord])
                removable = self.select(self.boundary_elements, all_tetra_ref[tetra_analize][vol_tetra == vol_tetra.min()][0])
                remov = rng.unite(remov,removable)

            self.create_vtk()
            import pdb; pdb.set_trace()
            self.remove_vol(removable)

    def find_bad_vol(self):
        index = 0
        point_type = np.concatenate((np.array([-1]), self.point_type))
        while not self.condition():
            #import pdb; pdb.set_trace()
            index = index + 1
            connect = self.connectivities(self.boundary_elements, 4)
            flag = point_type[connect]
            np.all(flag == 3,axis=1)



            import pdb; pdb.set_trace()

            all_tetra_ref = np.arange(connect.shape[0])
            num_tetra_bnodes = np.sum(np.isin(connect,np.array(self.boundary_nodes)), axis=1)
            four_bnodes_tetra = (num_tetra_bnodes == 4)
            three_bnodes_tetra = (num_tetra_bnodes == 3)
            # treating tetrahedron with 4 nodes on the boundary those tetrahedron
            # will be removed as no node are exclusive to these entities
            rm1 = all_tetra_ref[four_bnodes_tetra][~np.any(np.logical_and(connect[four_bnodes_tetra],flag[four_bnodes_tetra] == 0),axis=1)]
            #~np.any(np.logical_and(connect[four_bnodes_tetra],flag[four_bnodes_tetra] == 0),axis=1)
            #import pdb; pdb.set_trace()

            # treating tetrahedron with 3 nodes on the boundary those tetrahedron
            # will be removed as no node are exclusive to these entities
            rmp = all_tetra_ref[three_bnodes_tetra]
            # ~np.isin(connect[three_bnodes_tetra],np.array(self.boundary_nodes))
            tag = np.any(flag[three_bnodes_tetra] == 0,  axis=1)
            rm2 = rmp[~tag]
            tetra_tag = rmp[tag]
            rm3 = tetra_tag[np.all(~np.isin(connect[tetra_tag],np.array(self.boundary_nodes)) == (flag[tetra_tag] == 0), axis=1)]
            # flag[three_bnodes_tetra]
            # rm2 = all_tetra_ref[four_bnodes_tetra]
            remove_tetra = np.sort(np.hstack((rm1,rm2,rm3)))
            remove_flag = np.isin(all_tetra_ref, remove_tetra)

            self.remove_vol(self.select(self.boundary_elements,remove_flag))

            self.create_vtk()
            #import pdb; pdb.set_trace()

    def smooth(self):
        point_type = np.concatenate((np.array([-1]), self.point_type))
        connect = self.connectivities(self.boundary_elements, 4)
        all_tetra_ref = np.arange(connect.shape[0])
        flag = point_type[connect]
        self.remove_vol(self.select(self.boundary_elements,np.all(flag == 3,axis=1)))

    def condition(self):
        #desired_bnodes = np.array(self.nodes)[self.point_type == 0]
        print("desired nodes:", self.desired_bnodes)
        print("boundary nodes:", np.array(self.boundary_nodes))
        print("nodes on boundary", self.desired_bnodes[np.isin(self.desired_bnodes,np.array(self.boundary_nodes))])
        print("missing", self.desired_bnodes[~np.isin(self.desired_bnodes,np.array(self.boundary_nodes))])
        #print(~desired_bnodes[np.isin(desired_bnodes,np.array(self.boundary_nodes))])
        return np.all(np.isin(self.desired_bnodes,np.array(self.boundary_nodes)))

    def select(self, range, vec):
        return rng.Range(np.array(range)[vec])

    def find_boundary_nodes(self):
        return self.adjs(self.faces, 0)

    def find_boundary_elements(self):
        #self.mtu = topo_util.MeshTopoUtil(self.mb)
        rnb = rng.Range()
        it = ([self.mtu.get_bridge_adjacencies(el, 2,3) for el in self.faces])

        for el in it:
            rnb = rng.unite(rnb, el)
        return rnb


    def create_vtk(self):
        for index,el in enumerate(self.all_elements):
            self.mb.tag_set_data(self.tag_handle, el, index)
        for index,el in enumerate(self.nodes):
            mp = self.point_type[index].astype(int)
            # import pdb; pdb.set_trace()
            self.mb.tag_set_data(self.tag_node, el, mp)
        meshset = self.mb.create_meshset()
        self.mb.add_entity(meshset, self.nodes)
        self.mb.write_file("delaunay0.vtk", [meshset])
        self.mb.write_file("delaunay1.vtk")

    def remove_vol(self, element):
        tetra_remove = rng.intersect(self.all_elements, element)
        tetra_remaining = rng.subtract(self.all_elements, tetra_remove)
        #import pdb; pdb.set_trace()
        faces = rng.subtract(self.adjs(tetra_remove, 2), self.adjs(tetra_remaining,2))
        edges = rng.subtract(self.adjs(tetra_remove, 1), self.adjs(tetra_remaining,1))
        nodes = rng.subtract(self.adjs(tetra_remove, 0), self.adjs(tetra_remaining,0))
        self.mb.delete_entity(tetra_remove)
        self.mb.delete_entity(faces)
        self.mb.delete_entity(edges)
        self.mb.delete_entity(nodes)
        self.all_elements = self.mb.get_entities_by_dimension(0, 3)
        self.nodes = self.mb.get_entities_by_dimension(0, 0)
        self.faces = self.skinner(self.all_elements)
        self.boundary_nodes = self.find_boundary_nodes()
        self.boundary_elements = self.find_boundary_elements()

    def adjs(self, range,x):
        rnb = rng.Range()
        #import pdb; pdb.set_trace()
        it = ([self.mb.get_adjacencies(el, x) for el in range])
        for el in it:
            rnb = rng.unite(rnb, el)
        return rnb

    def connectivities(self,range, i):
        return self.mb.get_connectivity(range).reshape((-1, i))

    def external_nodes(self):
        range = self.faces
        return np.unique(np.array(self.mb.get_connectivity(range)))

    def skinner(self, range):
        skin = sk.Skinner(self.mb)
        return skin.find_skin(self.rs,range)

class DelaunayView(object):
    def __init__(self, coords, elements, neighbors, center):
        self.mb = core.Core()
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        skin = sk.Skinner(self.mb)
        verts = self.mb.create_vertices(coords)
        import pdb; pdb.set_trace()
        rs = self.mb.get_root_set()
        #import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()
        import pdb; pdb.set_trace()
        elements = elements + np.ones(elements.shape)
        tag_handle = self.mb.tag_get_handle("Teste", 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, create_if_missing = True)
        for el in elements:
            tetra = self.mb.create_element(types.MBTET, el.ravel().astype("uint64"))
        #self.mtu.construct_aentities(verts)
        elements = self.mb.get_entities_by_dimension(0, 3)
        for index,el in enumerate(elements):
            self.mb.tag_set_data(tag_handle,el,index)
        self.mb.write_file("delaunayVIEW.vtk")

class Delaunay(object):
    """Class for Visualization using moab of a Delaunay triangulation using SciPy."""

    def __init__(self, coords, elements, neighbors, center):
        #super(DelaunayView, self).__init__()
        self.mb = core.Core()
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        skin = sk.Skinner(self.mb)
        verts = self.mb.create_vertices(coords)
        rs = self.mb.get_root_set()
        #import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()
        elements = elements + np.ones(elements.shape)

        tag_handle = self.mb.tag_get_handle("Teste", 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, create_if_missing = True)


        for el in elements:
            tetra = self.mb.create_element(types.MBTET, el.ravel().astype("uint64"))


        #self.mtu.construct_aentities(verts)
        elements = self.mb.get_entities_by_dimension(0, 3)

        for index,el in enumerate(elements):
            self.mb.tag_set_data(tag_handle,el,index)
        #import pdb; pdb.set_trace()
        bfaces = skin.find_skin(0, elements)
        print(bfaces)
        adj = np.array([(self.mb.get_connectivity(bface)) for bface in bfaces])
        adj = adj.reshape((len(adj), 3)).astype(int) - np.ones(adj.shape)
        missing = np.setdiff1d(np.where(~np.isin(neighbors,adj)), center)
        import pdb; pdb.set_trace()
        # finding boundary tetrahedron
        self.mb.write_file("delaunay01.vtk")
        emp = rng.Range()
        boundary_tetrahedron = rng.Range()
        for el in elements:
            boundary_intersect = rng.intersect(bfaces, self.mb.get_adjacencies(el, 2))
            if boundary_intersect is not emp:
                boundary_tetrahedron = rng.unite(boundary_tetrahedron, rng.Range(el))


        # for el in elements:
        #     #import pdb; pdb.set_trace()
        #     local = rng.Range(el)
        #     boundary_intersect = rng.intersect(bfaces, self.mb.get_adjacencies(local, 2))
        #     if boundary_intersect is not emp:
        #         # import pdb; pdb.set_trace()
        #         face_con = self.mb.get_adjacencies(boundary_intersect,0)
        #         el_con = self.mb.get_adjacencies(local, 0)
        #             #import pdb; pdb.set_trace()
        #         inside_node = int(rng.subtract(el_con, face_con)[0] - 1)
        #
        #         #inside_node = int(np.setdiff1d(el_con, face_con)[0] - 1)
        #         #import pdb; pdb.set_trace()
        #         #check if inside node is missing
        #         is_missing = bool(np.isin(inside_node, missing))
        #         if is_missing is True:
        #             boundary_tetrahedron = rng.unite(boundary_tetrahedron, local)
                #import pdb; pdb.set_trace()


        for el in boundary_tetrahedron:
            self.mb.tag_set_data(tag_handle,el,2)
            print(self.mb.get_connectivity(el))
        #print(boundary_tetrahedron)
        import pdb; pdb.set_trace()

        #self.mb.delete_entity(boundary_tetrahedron)
        print(boundary_tetrahedron)
        self.mb.write_file("delaunay02.vtk")
