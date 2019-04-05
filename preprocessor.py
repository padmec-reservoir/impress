# Run preprocessor

import time
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
import preprocessor.geoUtil.geoTools as gtool

start = time.time()
M = msh("semi.msh", dim = 3)
end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
