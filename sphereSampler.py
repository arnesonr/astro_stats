import numpy as np
# from matplotlib.mlab import griddata
from scipy.interpolate import griddata

class sphereSampler(object):

    def __init__(self, radius = 1.0):
        """
        Class for creating and analysing samples on a sphere of radius r.

        :param radius (optional): radius of sphere, default: 1
        """
        self.r = radius

    def getSamples(self, n):
        """
        Get n uniform sample points from a sphere. Returns angles theta,
        phi of samples.

        :param n: number of samples
        :return theta, phi: angles of samples
        """
        phi = 2 * np.pi * np.random.rand(n)
        theta = np.arccos(2 * np.random.rand(n) - 1)
        return theta, phi

    def getCartesianCoords(self, theta, phi):
        """
        Routine to transform spherical into cartesian coordinates.

        :param theta: angles theta
        :param phi: angles phi
        :return x, y, z: cartesian coordinates of sample
        """
        x = self.r * (np.sin(theta) * np.cos(phi))
        y = self.r * (np.sin(theta) * np.sin(phi))
        z = self.r * np.cos(theta)
        return x, y, z

    def getSphericalCoords(self, x, y, z):
        """
        Routine to transform cartesian into spherical coordinates.

        :param x: x value
        :param y: y value
        :param z: z value
        :return r, theta, phi: spherical coordinates of samples
        """
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z/r)
        phi = np.arctan(y/x)
        return r, theta, phi

    def getF(self, theta, phi):
        """
        Evaluate dummy function F at sample points theta and phi.

        :param theta: angles theta
        :param phi: angles phi
        :return f: value of dummy function at (theta, phi)
        """
        return np.absolute(theta-.5*np.pi)
        # return np.cos(theta)+np.sin(theta)

    def getInterpolatedGrid(self, x, y, f, nbins = 21, method = 'linear'):
        """
        Put random samples on regular cartesian grid and interpolate the
        values of f.

        :param x: x coordinates of samples
        :param y: y coordinates of samples
        :param f: function values of samples
        :param nbins: number of bins per axis
        :param method: interpolation method of griddata
        :return X, Y, F: x and y coordinates of grid, interpolated function
        values
        """
        xx = np.linspace(-self.r, self.r, nbins)
        yy = np.linspace(-self.r, self.r, nbins)
        X, Y = np.meshgrid(xx, yy)
        grid = griddata(points = np.vstack((x, y)).T, values = f, xi = (X, Y), method=method)
        return X, Y, grid
