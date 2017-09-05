import numpy as n
from matplotlib import pyplot as p
import scipy.stats as stats
from astropy.io import fits
import scipy.integrate as integ
import scipy.optimize as opt

def gaussianPDF(mu, sigma, x):
    z=(x - mu) / (sigma)
    return 1./(sigma * n.sqrt(2 * n.pi)) * n.exp( - z**2/2 )

def mll(theta, data):
    mu = theta[0]
    sigma = theta[1]
    pxs = n.log(gaussianPDF(mu, sigma, data))
    return -n.sum((pxs))

def main():
	hdulist = fits.open('ABELL_3266.fits')
	#data in form [objID,run,rerun,camcol,field,obj,type,ra,dec,u,g,r,i,z,
	#err_u,err_g,err_r,err_i,err_z]
	data = hdulist[1].data
	hdr = hdulist[1].header
	hdulist.close()
	#a) plot histogram of g-r
	vel = data['Vel']
	e_vel = data.field('e_Vel')
	#make a guess on the average and sigma
	initial=[10000.,500.]
	#optimize the fit
	v_fit = opt.minimize(mll,initial,args=(vel,))

	print 'Fitted mu:', v_fit.x[0]
	print 'Fitted sigma:', v_fit.x[1]
 	#we measured the total error which factors in the errors of each velocity
 	#i.e. sigma_i^2 = sigma_o^2 + e_i^2
 	#calculate what sigma_i is
	sigma_o = res.x[1]
	sigma_i = n.sqrt(sigma_o**2 + n.sum(e_vel**2))

	e_vel_fit = opt.minimize(mll,[500.,20.],args=(e_vel,))

if __name__ == '__main__':
    main()