from pymoab import core, types, rng, topo_util
from pymoab import skinner as sk
import numpy as np


def volume(p1, p2, p3, p4):
    d1 = p2 - p1
    d2 = p3 - p1
    d3 = p4 - p1
    matrix = np.vstack((d1, d2, d3))
    return (1/6)*np.linalg.det(matrix)



class DelaunaySingle(object):
    def __init__(self, coords, elements, point_type, center):
        self.mb = core.Core()
        self.point_type = point_type
        verts = self.mb.create_vertices(coords)
        self.rs = self.mb.get_root_set()
        #import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()
        elements = elements + np.ones(elements.shape)
        self.tag_handle = self.mb.tag_get_handle("Teste", 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, create_if_missing = True)
        self.tag_node = self.mb.tag_get_handle("Nodes", 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, create_if_missing = True)

        for el in elements:
            tetra = self.mb.create_element(types.MBTET, el.ravel().astype("uint64"))

        # self.mtu = topo_util.MeshTopoUtil(self.mb)
        # self.mtu.construct_aentities(verts)
        #self.mtu.construct_aentities(verts)
        #self.elements = self.mb.get_entities_by_dimension(0, 3)
        self.all_elements = self.mb.get_entities_by_dimension(0, 3)
        self.nodes = self.mb.get_entities_by_dimension(0, 0)
        self.faces = self.skinner(self.all_elements)
        self.boundary_nodes = self.find_boundary_nodes()
        self.boundary_elements = self.find_boundary_tetra()
        self.remove_vol(self.all_elements[0:120])
        self.create_vtk()


    def find_boundary_nodes(self):
        return self.adjs(self.faces, 0)

    def find_boundary_elements(self):
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        import pdb; pdb.set_trace()
        self.mtu.bridge_adjacencies()
        return self.adjs(self.faces, 0)


    def create_vtk(self):
        for index,el in enumerate(self.all_elements):
            self.mb.tag_set_data(self.tag_handle, el, index)
        for index,el in enumerate(self.nodes):
            self.mb.tag_set_data(self.tag_node, el, self.point_type[index].astype(int))
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

    def adjs(self, range,x):
        rnb = rng.Range()
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
