
# Run preprocessor

# docker run -t -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"
# docker run -t -it -v  /home/arturcastiel/projetos:/pytest gabrielmmats/impress bash -c "cd /pytest; bash"
# docker run -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"
# %load_ext autoreload
# %autoreload 2

from pymoab import core, types, rng, topo_util, skinner, tag
import numpy as np
import time
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
from preprocessor.meshHandle.finescaleMesh import FineScaleMesh
# from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual
import preprocessor.geoUtil.geoTools as gtool

start = time.time()
M = msh('mesh/eyem.msh', dim=2)

faces = M.core.all_faces[:]
ifaces = M.coarse.all_interface_faces[:]
iface_tag = M.core.mb.tag_get_handle("iface", size=1, tag_type=types.MB_TYPE_INTEGER, 
                                    storage_type=types.MB_TAG_SPARSE, create_if_missing=True)
M.core.mb.tag_set_data(iface_tag, faces[ifaces], np.ones(ifaces.size, dtype=int))

M.core.print(file='preprocessed_eyem', extension='.vtk')

end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
