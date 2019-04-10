# Procedure to get the relevant vertices for the geometric module in each type of elements

from pymoab import core, topo_util, rng


mb = core.Core()
mesh_topo_util = topo_util.MeshTopoUtil(mb)

# Piramide


def get_pyramid_adjacencies(vertices, edges):

    # Creating used ranges
    aux1 = rng.Range()
    aux2 = rng.Range()
    base_node = rng.Range()
    top_node = rng.Range()

    # Discovering the top node of the pyramid
    for e in vertex_handles:
        adj = mesh_topo_util.get_bridge_adjacencies(e, 1, 0)
        if len(adj) == 4:
            top_node.insert(e)

    # Defines the base nodes
    base_nodes = rng.subtract(vertex_handles, top_node)

    # Getting the base adjacencies of the point of reference
    aux = mesh_topo_util.get_bridge_adjacencies(base_nodes[0], 1, 0)
    aux = rng.subtract(top_node, aux)

    aux1.insert(aux[0])
    aux2.insert(aux[1])
    base_node.insert(base_nodes[0])
    coord1 = mbcore.get_coords(base_node)
    coord2 = mbcore.get_coords(aux1)
    coord3 = mbcore.get_coords(aux2)
    coord4 = mbcore.get_coords(top_node)
    order = ([coord1, coord2, coord3, coord4])

    return order

def get_hexahedro_adjacencies(vertices, edges):

    # Creating used ranges
    aux1 = rng.Range()
    aux2 = rng.Range()
    base_nodes = rng.Range()
    top_node = rng.Range()

    # Discovering the top node of the pyramid
    for e in vertex_handles:
        adj = mesh_topo_util.get_bridge_adjacencies(e, 1, 0)
        if len(adj) == 4:
            top_node.insert(e)

    # Defines the base nodes
    base_nodes = rng.subtract(vertex_handles, top_node)

    # Getting the base adjacencies of the point of reference
    aux = mesh_topo_util.get_bridge_adjacencies(base_nodes[0], 1, 0)
    aux = rng.subtract(top_node, aux)

    aux1.insert(aux[0])
    aux2.insert(aux[1])
    base_node.insert(base_nodes[0])
    coord1 = mbcore.get_coords(base_node)
    coord2 = mbcore.get_coords(aux1)
    coord3 = mbcore.get_coords(aux2)
    coord4 = mbcore.get_coords(top_node)
    order = ([coord1, coord2, coord3, coord4])

    return order
