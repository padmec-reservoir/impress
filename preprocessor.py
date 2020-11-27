
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
# from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual
import preprocessor.geoUtil.geoTools as gtool

start = time.time()
M = msh('mesh/fivespot105u.msh', dim=2)
# M = msh('mesh/icecream.msh', dim=3)
M.core.print(file='preprocessed_icream', extension='.vtk')
for n, coarse_elem in zip(np.arange(len(M.coarse.elements)), M.coarse.elements):
    print("writing element {}".format(n+1))
    coarse_elem.core.print(file='preprocessed_icecream_{}'.format(n), extension='.vtk')
# M.coarse.elements[0].core.print(file='teste_coarse_element', extension='.vtk')
end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
