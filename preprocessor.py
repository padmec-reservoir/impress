
# Run preprocessor

# docker run -t -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"
# docker run -t -it -v  /home/arturcastiel/projetos:/pytest gabrielmmats/impress bash -c "cd /pytest; bash"
# docker run -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"
# %load_ext autoreload
# %autoreload 2

from pymoab import core, types, rng, topo_util, skinner, tag
import numpy as np
import time
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS
from preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh
# from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual
import preprocessor.geoUtil.geoTools as gtool
from preprocessor.element_order.MPFAD2DOrdering import MPFAD2DOrdering

start = time.time()
M = msh('mesh/MeshSkewed_mod_02_2x2.msh', dim=2)

some_edges = M.nodes.bridge_adjacencies(8, "nodes", "edges")
some_edges_ord = M.nodes.bridge_adjacencies(8, "nodes", "edges", ordering_inst=MPFAD2DOrdering(M.edges, "faces"))

print("Edges: {}".format(some_edges))
print("Edges after sort: {}".format(some_edges_ord))

# edges_set = M.core.mb.create_meshset()
# M.core.mb.add_entities(edges_set, M.core.mb.get_entities_by_dimension(0, 1))
# M.core.mb.write_file("MeshSkewed_mod_02_2x2_edges.vtk", output_sets=[edges_set])

# order = MPFAD2DOrdering(M.nodes, "faces")
# ord_nodes = M.nodes.bridge_adjacencies(M.nodes.all[2], "edges", "nodes", ordering_inst=order)
# print(ord_nodes)

end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
