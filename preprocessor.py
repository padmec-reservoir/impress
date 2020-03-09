
# Run preprocessor

# docker run -t -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"
# docker run -t -it -v  /home/arturcastiel/projetos:/pytest gabrielmmats/impress bash -c "cd /pytest; bash"
# docker run -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"
# %load_ext autoreload
# %autoreload 2
import numpy as np
import time
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual
import preprocessor.geoUtil.geoTools as gtool

start = time.time()
M = msh('mesh/40.h5m', dim = 3)
# M = msh('mesh/malha03.msh', dim = 2)
end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
#l = dual(M)
