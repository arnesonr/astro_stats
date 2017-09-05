import numpy as n
from matplotlib import pyplot as p
import scipy.stats as stats
from astropy.io import fits
import scipy.integrate as integ
import scipy.optimize as opt
from astropy.io import ascii

x    = n.array([1, 2, 3, 4, 5])
y    = n.array([1, 2, 3, 6, 4])
yerr = n.array([[0.1, 0.1, 0.1, 0.1, 0.1],
                 [1, 1, 1, 1, 1],
                 [1, 1, 1, .2, 2],
                 [1, 1, 1, 2, .2]])

pol, cov = n.polyfit(x, y, 1, w=1/yerr[1], cov=True)  # coefficients and covariance matrix
yfit = n.polyval(pol, x)          # evaluate the polynomial at x

perr = n.sqrt(n.diag(cov))     # standard-deviation estimates for each coefficient
R2 = n.corrcoef(x, y)[0, 1]**2  # coefficient of determination between x and y
resid = y - yfit
chi2red = n.sum((resid/yerr)**2)/(y.size - 2)
p.figure(figsize=(10, 5))
p.errorbar(x, y, yerr=yerr, fmt = 'bo', ecolor='b', capsize=0)
p.plot(x, yfit, 'r', linewidth=3, color=[1, 0, 0, .5])
p.xlabel('x', fontsize=16)
p.ylabel('y', fontsize=16)
p.show()

data = ascii.read('spectrum.dat')
wave=data['Wavelength']
flux=data['Flux']
p.figure(figsize=(10, 5))
p.step(wave,flux,alpha=0.3,where='mid',lw=3,label='sdss spectrum')
p.xlabel('wavelength', fontsize=16)
p.ylabel('Flux', fontsize=16)
p.xlim(6780, 6860)
p.ylim(8, 14)
p.show()

def lorentzian(x, A1, A2, p01, p02, ):
    A1, p01, w, offset = p
    u1 = (p01 - x) / w1
    A2, p02, w, offset = p
    u2 = (p02 - x) / w2
    return (A1 * 1.0/(1.0+u1**2) + offset) + (A2 * 1.0/(1.0+u2**2) + offset) 

