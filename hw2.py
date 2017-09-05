from astropy.io import fits
from matplotlib import pyplot as p
import numpy as n
import scipy.integrate as integ
import scipy.stats as stats
from numpy import *
import scipy.optimize as opt

def bootstrap(x,y):
	"""
  	PURPOSE: draws a bootstrap sample from x and y arrays

  	ARGUMENTS:
    	x: 1d array
    	y: (matching) 1d array
  	RETURNS: bootstrapped arrays x,y
  	"""
  	num = n.size(x)
  	boot_x = n.zeros(num)
  	boot_y = n.zeros(num)
  	for i in range(num):
  		j = n.random.randint(num)
  		boot_x[i] = x[j]
  		boot_y[i] = y[j]
	return [boot_x, boot_y]

def main():
	#
	#1.
	#
	hdulist = fits.open('test_correlation.fits')
	#data in form 
	data = hdulist[0].data
	hdulist.close()
	pairs = n.shape(data)[0]
	#a)
	pairs = n.shape(data)[0]
	#seperate all the x's and y's
	x = data[:,0]
	y = data[:,1]
	#calculate the Pearson's correlation coefficient and p value
	r_p, p_p = stats.pearsonr(x,y)
	#compute the Spearman's correlation coefficient and p value
	r_s, p_s = stats.spearmanr(x,y)
	#compute the Kendall tau coefficient and p value
	tau, p_t = stats.kendalltau(x,y)
	#the significance levels are:
	print "Significance (Pearson)", p_p 
	print "Significane (Spearman)", p_s 
	print "Significane (Kendall):", p_t
	#the significance using the Spearman's r coefficient is:
	#t=(r_s*n.sqrt(pairs-2))/(n.sqrt(1-r_s**2))
	#b)
	b1 = bootstrap(x,y)
	p.plot(x,y,'k.') #original sample
	p.plot(b1[0],b1[1], 'rx') #bootstrap sample
	#c)
	boot_mc = n.zeros((1000,2,200))
	for i in range(1000):
		boot_mc[i] = bootstrap(x,y)
	r_p = n.zeros(1000)
	r_s = n.zeros(1000)
	tau = n.zeros(1000)
	for i in range(1000):
		r_p[i] = stats.pearsonr(boot_mc[i][0],boot_mc[i][1])[0]
		r_s[i] = stats.spearmanr(boot_mc[i][0],boot_mc[i][1])[0]
		tau[i] = stats.kendalltau(boot_mc[i][0],boot_mc[i][1])[0]
	#d)
	#pearsons coefficient
	p.hist(r_p,bins=100)
	#expectation value of r is rho*(1-(1-rho^2)/N)
	#from the covariance matrix rho = 0.9
	Er = 0.9*(1.-((1.-.9**2)/n.size(x)))
	p.plot((Er,Er),(0,40),"k--",label='Expectation Value')
	p.title('Pearson Coefficient')
	p.legend(loc=2)
	#Spearmans coeffiient
	p.hist(r_s,bins=100)
	p.plot((Er,Er),(0,40),"k--",label='Expectation Value')
	p.title('Pearson Coefficient')
	p.legend(loc=2)
	#Kendall
	#the expectation value of tau is 2/pi*arcsin(rho)
	Etau = (2./n.pi)*n.arcsin(0.9)
	p.hist(tau,bins=100)
	p.plot((Etau,Etau),(0,35),'k--')

	#e)
	#add in some outliers
	x2 = n.concatenate((x,[-0.6,-0.76,-2.0]))
	y2 = n.concatenate((y,[-2.1,0.92,-5.0]))
	boot_mc2 = n.zeros((1000,2,203))
	for i in range(1000):
		boot_mc2[i] = bootstrap(x2,y2)
	r_p2 = n.zeros(1000)
	r_s2 = n.zeros(1000)
	tau2 = n.zeros(1000)
	for i in range(1000):
		r_p2[i] = stats.pearsonr(boot_mc2[i][0],boot_mc2[i][1])[0]
		r_s2[i] = stats.spearmanr(boot_mc2[i][0],boot_mc2[i][1])[0]
		tau2[i] = stats.kendalltau(boot_mc2[i][0],boot_mc2[i][1])[0]
	p.hist(r_p2,bins=100)
	p.plot((.9,.9),(0,35),"k--")
	p.hist(r_s2,bins=100)
	p.plot((.9,.9),(0,35),"k--")
	p.hist(r_p,bins=100)
	p.plot((.9,.9),(0,35),"k--")
	#it lowers the correlation coefficients (makes them less correlated)


	#2.
	#a)
	ew_data =loadtxt('EW_Ouchi2007.dat')
	p.hist(ew_data,bins=50,normed=1)

	#b)paper
	#c)
	def Lya_PDF(Wo, EW):
		return 1./(Wo*n.exp(-50./Wo)) * n.exp( - EW/Wo)

	def logl(Wo,data):
		pxs = n.log(Lya_PDF(Wo, data))
		return n.sum((pxs))

	def mll(Wo,data):
		pxs = n.log(Lya_PDF(Wo, data))
		return -n.sum((pxs))

	#use ML to find best estimate for Wo
	Wos = n.arange(1.0,350.,0.01)
	lls = n.zeros(n.size(Wos))
	for i in range(n.size(Wos)):
			lls[i] = logl(Wos[i],ew_data)
	Wo = Wos[n.argmax(lls)]
	p.plot(Wos,lls,'k-')
	p.xlabel('Wo')
	p.ylabel('Ln(Likelihood)')
	#d)
	#graphic method we want to find where Ln(L) is lower by 1/2
	def find_nearest_element(array,value,index=False):
		"""Return closest value in array, or its index"""
		idx = n.abs(array-value).argmin()
		return (idx,array.flat[idx]) if index else array.flat[idx]
	max_index = find_nearest_element(lls,(n.max(lls)),index=True)[0]
	#Wo minus 1 sigma
	Wo_msigma = Wos[find_nearest_element(lls[0:max_index],(n.max(lls)-.5),index=True)[0]]
	Wo_psigma = Wos[find_nearest_element(lls[max_index:5000],(n.max(lls)-.5),index=True)[0]+max_index]
		#sigma of Wo
	p.plot(Wos[0:3500],lls[0:3500],'k-')
	p.plot((Wo_msigma,Wo_msigma),(-180,-182),'r--')
	p.plot((Wo_psigma,Wo_psigma),(-180,-182),'r--')
	p.plot((Wo,Wo),(-180,-182),'b.-.')
	#using the bootstrap method
	def bootstrap_wo(data):
		"""
		PURPOSE: runs a ML bootstrap sample to estimate error on Wo

		ARGUMENTS:
			data: data to bootstrap			

		RETURNS: ML estimates of Wo
		"""
		ew_boot = n.zeros(n.size(data))
		for i in range(n.size(data)):
			j = n.random.randint(n.size(data))
			ew_boot[i] = data[j]
		return ew_boot
	#need n(log(n)**2) ~ 445 bootstraps
	#analytical result for Wo is n.mean(ew_data) - 50.0
	Wo_ML = n.zeros(500)
	for k in range(500):
		bootstrap_data = bootstrap_wo(ew_data)
		Wo_ML[k] = (n.sum(bootstrap_data)/n.size(bootstrap_data)-50.0
	p.hist(Wo_ML,bins=100)
	print "Wo_ML Std dev:", n.std(Wo_ML)
if __name__ == '__main__':
    main()