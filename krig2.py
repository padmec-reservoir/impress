from pykrige.uk import UniversalKriging
from pykrige.ok import OrdinaryKriging
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


data = np.array([[0.3, 1.2, 0.47],
                 [1.9, 0.6, 0.56],
                 [1.1, 3.2, 0.74],
                 [3.3, 4.4, 1.47],
                 [2.3, 5.4, 1.67],
                 [2,4,1.7],
                 [4.7, 3.8, 1.74]])

nx = 3
ny = 3
gridx = np.arange(0.0, 5.5, 1/nx)
gridy = np.arange(0.0, 5.5, 1/ny)

# Create the ordinary kriging object. Required inputs are the X-coordinates of
# the data points, the Y-coordinates of the data points, and the Z-values of the
# data points. Variogram is handled as in the ordinary kriging case.
# drift_terms is a list of the drift terms to include; currently supported terms
# are 'regional_linear', 'point_log', and 'external_Z'. Refer to
# UniversalKriging.__doc__ for more information.
UK = UniversalKriging(data[:, 0], data[:, 1], data[:, 2], variogram_model='linear',
                      drift_terms=['regional_linear'])

OK = OrdinaryKriging(data[:, 0], data[:, 1], data[:, 2], variogram_model='linear',
                     verbose=False, enable_plotting=False)

zordinary, ss = OK.execute('grid', gridx, gridy)
# Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
# grid of points, on a masked rectangular grid of points, or with arbitrary points.
# (See UniversalKriging.__doc__ for more information.)
z, ss = UK.execute('grid', gridx, gridy)


xx, yy = np.meshgrid(gridx, gridy)

# dzz = z.ravel()
# dzo = zordinary.ravel()

dzo = zordinary.ravel()
dzz = z.ravel()
dxx = xx.ravel()
dyy = yy.ravel()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
#ax.plot_surface(dxx, yy, dzz.data, color='b')


#

# ax.scatter(gridx, gridy, z, c='b', marker='x')
#ax.scatter(dxx, dyy, dzz, c='b', marker='o')
# ax.scatter(gridx, gridy, z, c='b', marker='x')
ax.scatter(dxx, dyy, dzo, c='b', marker='x')
ax.scatter(data[:,0], data[:,1], data[:,2], c='r', marker='o')
#
# ax.scatter(dxx, dyy, dzz, c='r', marker='o')

#
plt.show()
