from pykrige.uk import UniversalKriging
from pykrige.ok3d import OrdinaryKriging3D
from pykrige.uk3d import UniversalKriging3D

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


data = np.array([[0.1, 0.1, 0.3, 0.9],
                                 [0.2, 0.1, 0.4, 0.8],
                                 [0.1, 0.3, 0.1, 0.9],
                                 [0.5, 0.4, 0.4, 0.5],
                                 [0.3, 0.3, 0.2, 0.7]])

gridx = np.arange(0.0, 0.6, 0.05)
gridy = np.arange(0.0, 0.6, 0.01)
gridz = np.arange(0.0, 0.6, 0.1)

# Create the ordinary kriging object. Required inputs are the X-coordinates of
# the data points, the Y-coordinates of the data points, and the Z-values of the
# data points. Variogram is handled as in the ordinary kriging case.
# drift_terms is a list of the drift terms to include; currently supported terms
# are 'regional_linear', 'point_log', and 'external_Z'. Refer to
# UniversalKriging.__doc__ for more information.
UK = UniversalKriging(data[:, 0], data[:, 1], data[:, 2], variogram_model='linear',
                      drift_terms=['regional_linear'])

# Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
# grid of points, on a masked rectangular grid of points, or with arbitrary points.
# (See UniversalKriging.__doc__ for more information.)
z, ss = UK.execute('points', gridx, gridy)


xx, yy = np.meshgrid(gridx, gridy)

dzz = z.ravel()
dxx = xx.ravel()
dyy = yy.ravel()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
#ax.plot_surface(dxx, yy, dzz.data, color='b')


#

ax.scatter(gridx, gridy, z, c='b', marker='x')
# ax.scatter(dxx, dyy, dzz, c='r', marker='o')
ax.scatter(data[:,0], data[:,1], data[:,2], c='r', marker='x')
#
# ax.scatter(dxx, dyy, dzz, c='r', marker='o')

#
plt.show()
