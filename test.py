from pymoab import core, types, rng, topo_util
from pymoab import skinner as sk

mb = core.Core()
mtu = topo_util.MeshTopoUtil(mb)
mb.load_file("40.h5m")
all_nodes = mb.get_entities_by_dimension(0, 0)
mtu.construct_aentities(all_nodes)
all_volumes = mb.get_entities_by_dimension(0, 3)
my_volumes = all_volumes[10000:20000]

my_tag = tag_get_handle("my_tag", 1, pymoab.types.MB_TYPE_INTEGER, pymoab.types.MB_TAG_DENSE, True, 0)
mb.tag_set_data(my_tag, all_volumes, 50)

my_meshset = mb.create_meshset()
mb.add_entities(my_meshset, my_volumes)

#print meshset
mb.write_file("test_file.h5m", [0])
