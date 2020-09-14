import numpy as np
M.dreams[:] = 0
M.dreams[np.unique(np.concatenate(l.coarse_edges))] = 3
