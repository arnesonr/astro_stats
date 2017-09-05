from matplotlib import pyplot as pl
import numpy as np

def plot_scatter(x, y, c = None, xlabel = '', ylabel = '', clabel = '', radius = None):
    """
    Scatter plot with optional color coding when using c.

    :param x: x values of samples
    :param y: y values of samples
    :param c (optional): function values of samples to be color coded,
    default: None
    :param xlabel (optional): label of x-axis, default: no label
    :param ylabel (optional): label of y-axis, default: no label
    :param radius (optional): overplot circle of given radius, default: None
    """
    if c != None:
        pl.scatter(x, y, c = c)
        pl.colorbar(label = clabel)
    else:
        pl.scatter(x, y)
    if xlabel != None:
        pl.xlabel(xlabel)
    pl.ylabel(ylabel)
    if radius != None:
        fig = pl.gcf()
        circle = pl.Circle((0,0),radius,fill=False)
        fig.gca().add_artist(circle)
        fig.gca().set_aspect(1.0)
    pl.show()

def plot_histogram(x, xlabel = '', ylabel = '', nbins = 10, nsamples = None):
    """
    Plot a histogram.

    :param x: x values of samples
    :param xlabel (optional): label of x-axis, default: no label
    :param ylabel (optional): label of y-axis, default: no label
    :param nbins (optional: number of bins for histogram, default: 10
    :param nsamples (optional): plot line of expected number of samples per
    bin, default: None
    """
    pl.hist(x, nbins)
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)
    if nsamples != None:
        pl.axhline(float(nsamples)/nbins, c = 'k')
    pl.show()

def plot_datagrid(X, Y, Z, xlabel = '', ylabel = ''):
    """

    :param X:
    :param Y:
    :param Z:
    :return:
    """
    pl.imshow(Z.T, extent=(np.amin(X),np.amax(X),np.amin(Y),np.amax(Y)), origin='lower', interpolation='none')
    pl.colorbar()
    pl.show()
