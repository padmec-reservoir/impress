# Run preprocessor

import time
from mspreprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
import mspreprocessor.geoUtil.geoTools as gtool

start = time.time()
M = msh("25.h5m", dim = 3)
end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
