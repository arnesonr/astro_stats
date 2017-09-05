import numpy as n
from matplotlib import pyplot as p
import scipy.stats as stats
from astropy.io import fits
import numpy.random as r
import scipy.optimize as opt
import sklearn.mixture as GMM

def pdf_model(x, theta):
    mu1, sig1, mu2, sig2, al_1 = theta
    return al_1*normpdf(x, mu1, sig1) + (1-al_1)*normpdf(x, mu2, sig2)

def log_likelihood_two_1d_gauss(theta, sample):
    	return -n.log(pdf_model(sample, theta)).sum()

def main():
	hdulist = fits.open('GC_NGC4365.fits')
	#hdulist.info()
	#data in form [objID,run,rerun,camcol,field,obj,type,ra,dec,u,g,r,i,z,
	#err_u,err_g,err_r,err_i,err_z]
	data = hdulist[1].data
	hdr = hdulist[1].header
	#hdr.keys()
	#data.columns.info()
	hdulist.close()
	ra = data['RAJ2000']
	dec = data['DEJ2000']
	gmag = data['gmag']
	imag = data['imag']
	gi = gmag - imag
	p.hist(gi,bins=50)
	
    # Initial guess of parameters [mu1, sig1, mu2, sig2, a_1]
    #by eye there looks like at least two gaussians
	theta0 = n.array([0.8,0.1,1.0,0.2,0.5])
	#use the minimize package of scipy optimize
	#there are a plethora of method options ‘Nelder-Mead’,‘Powell’,‘CG’,
	#‘BFGS’,‘Newton-CG’,‘Anneal',‘L-BFGS-B’,‘TNC’,‘COBYLA’,‘SLSQP’,‘dogleg’,
	#‘trust-ncg’
	#res = opt.minimize(log_likelihood_two_1d_gauss, x0=theta0, args=(gi,), method='BFGS')
	#test different amount of components
	means = list()
	weights = list()
	covars = list()
	aics = list()
	for i in range(5):
		model = GMM(i+1)
		model.fit(gi)
		means.append(model.means_.flatten())
		#weights
		weights.append(model.weights_.flatten())
		#covariances
		covars.append(model.covars_.flatten())
		M_best = model
		aics.append(M_best.aic(gi))
	#the GCs are best fit with just 2 components

	model = GMM(2)
	model.fit(gi)
	#print the means
	model.means_.flatten()
	#weights
	model.weights_.flatten()
	#covariances
	model.covars_.flatten()
	M_best = model
	x = n.arange(.6, 1.3,.01)
	logprob,respons = M_best.score_samples(x)
	pdf = n.exp(logprob)
	pdf_individual = respons * pdf[:,n.newaxis]
	#plot the figure
	fig, ax = p.subplots(1,1,figsize=(10,7))
	ax.hist(gi,50,normed=True,histtype='stepfilled',alpha=0.4,label='data')
	ax.plot(x,pdf,'-k',lw=3,label="Best fit")
	ax.plot(x,pdf_individual[:,0],'--',c='k',label='component 1')
	ax.plot(x,pdf_individual[:,1],'-',c='k',label='component 2')
	ax.legend()
	ax.set_xlable('x',fontsize=20)
	#calculate the Akaike information criteria for the model
	M_best.aic(gi)
	
	#get the component each object likely belongs to	
	component = M_best.predict(gi)
	ra1 = ra[n.where(component==0)]
	ra2 = ra[n.where(component==1)]
	dec1 = dec[n.where(component==0)]
	dec2 = dec[n.where(component==1)]
	p.plot(ra1,dec1,'b.',label='component 1')
	p.plot(ra2,dec2,'r.',label='component 2')
	p.hist(ra1,bins=50,alpha=.4,color='red')
	p.hist(ra2, bins=50,alpha=.4,color='blue')
	p.hist(dec1,bins=50,alpha=.4,color='red')
	p.hist(dec2, bins=50,alpha=.4,color='blue')
	#the GC are distributed differently

if __name__ == '__main__':
    main()