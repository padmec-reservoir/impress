
# Run preprocessor

# docker run -t -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"

# docker run -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"
# %load_ext autoreload
# %autoreload 2
import time
from preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh
from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual
#from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual2

#import preprocessor.geoUtil.geoTools as gtool
#import sys
#import imp



#print(sys.path)
#sys.path.append('/mesh')
#foobar = imp.load_source('20.h5m', '/mesh')

start = time.time()
M = msh('mesh/40.h5m', dim=3)
# M = msh('mesh/malha03.msh', dim = 2)


end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
#l = dual(M)


#start = time.time(); l = dual(M); end = time.time(); print("The preprocessing step lasted {0}s".format(end-start))
#end = time.time()

# start = time.time()
# l = dual2(M)
# print("The preprocessing step lasted {0}s".format(end-start))
# end = time.time()
