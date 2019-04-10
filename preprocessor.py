# Run preprocessor

import time
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
import preprocessor.geoUtil.geoTools as gtool
#import sys
#import imp

#print(sys.path)
#sys.path.append('/mesh')
#foobar = imp.load_source('20.h5m', '/mesh')

start = time.time()
M = msh('20.h5m', dim = 3)
end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
