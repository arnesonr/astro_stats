import numpy as np

def sample_sphere(n):
    """
    Routine to uniformly sample n points from a sphere. Returns angles theta,
    phi of samples.

    :param n: number of samples
    :return theta, phi: angles of samples
    """
    theta = 2 * np.pi * np.random.rand(n)
    phi = np.arccos(2 * np.random.rand(n) - 1)
    return theta, phi

def get_cartesian_coords(phi, theta, r = 1):
    """
    Routine to transform spherical into cartesian coordinates.

    :param phi: angles phi
    :param theta: angles theta
    :param r (optional): radius of sphere, default: 1
    :return x, y, z: cartesian coordinates of sample
    """
    x = r * (np.sin(phi) * np.cos(theta))
    y = r * (np.sin(phi) * np.sin(theta))
    z = r * np.cos(phi)
    return x, y, z
