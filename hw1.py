from astropy.io import fits
from matplotlib import pyplot as p
import numpy as n
import scipy.integrate as integ
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.stats as stats

def median(array):
	"""
  	PURPOSE: Calculate the median value of an 1d array

  	ARGUMENTS:
    	array: 1d array

  	RETURNS: median value
  	"""
  	#sort the array from low to high
	array.sort()
	#find the middle of the array
	mid = n.size(array)/2
	#if there is an even number of elements return the average of the middle two
	if n.size(array) % 2 == 0:
		return (array[mid-1] + array[mid]) / 2.0
	#otherwise return the middle of the array
	else:
		return array[mid]

def PDF(p,var):
	"""
  	PURPOSE: Calculate a value of x taken from p

  	ARGUMENTS:
    	p: function, A Python function or method to integrate
    	var: 1d increasing array, the range of the random variable x

  	RETURNS: a random number x generated according to the pdf
  	"""
  	u = n.random.uniform()
	#solve F(x) = u where F is the cdf of the given pdf
	#F = int[-inf,x] (pdf)
	#initialize F
	F = 0.0
	delta_x = 1.0 #initial increment to increase x
	x = n.min(var) #initial value to return
	while(F < u):
		Fl = F #the lower value of F
		F,Ferr = integ.quad(p, n.min(var),x)
		x += delta_x
	#refine x value i times depending on precision wanted
	for i in range(1,6):
		F = Fl
		x -= delta_x
		delta_x = 1.0/(10.**i)
		while(F < u):
			Fl = F
			F,Ferr = integ.quad(p, n.min(var),x)
			x += delta_x
	return x

def PDF2(p,var):
	"""
	My 2nd, more efficient version of randomly drawing from a pdf
	
  	PURPOSE: Calculate a value of x taken from p

  	ARGUMENTS:
    	p: function, A Python function or method to integrate
    	var: 1d increasing array, the range of the random variable var

  	RETURNS: a random number x generated according to the pdf
  	"""
  	switch = 0
  	while switch == 0:
	  	#draw two random numbers between the range of the pdf
	  	ux = n.random.uniform(n.min(var),n.max(var))
	  	uy = n.random.uniform(n.min(p(var)),n.max(p(var)))
	  	u = n.random.uniform()
		#decide if the random number, ux, should be kept
		#keep it if uy is less than p(ux)
		if uy <= p(ux):
			switch = 1
			return ux

def MC(k,m,p,var):
	"""
  	PURPOSE: Runs a Monte Carlo sampling of the pdf p

  	ARGUMENTS:
  		k: number of samples
  		m: data points in each sample
    	p: function, A Python function or method to integrate
    	var: 1d increasing array, the range of the random variable x

  	RETURNS: the matrix of random numbers, the sample means and stddev
  	"""
	d_k = n.zeros((k,m))
	for i in range(0,k):
		for j in range(0,m):
			d_k[i][j] = PDF2(p,var)
	d_k_means = n.zeros(k)
	d_k_std = n.zeros(k)
	for i in range(0,k):
		d_k_means[i] = n.mean(d_k[i][:])
		d_k_std[i] = n.std(d_k[i][:])
	return d_k, d_k_means, d_k_std

def main():
	#
	#1.
	#
	hdulist = fits.open('sdss.fits')
	#data in form [objID,run,rerun,camcol,field,obj,type,ra,dec,u,g,r,i,z,
	#err_u,err_g,err_r,err_i,err_z]
	data = hdulist[1].data
	hdr = hdulist[1].header
	hdulist.close()
	#a) plot histogram of g-r
	g = data.field('g')
	r = data.field('r')
	gmr = g-r
	c, bins, patches = p.hist(gmr,bins=30,log=True)
	p.xlabel('g-r')
	p.ylabel('number of objects')
	p.title('g-r Histogram')
	p.savefig('g_r_hist')
	p.close()
	#b)
	#median function is defined above
	#I did not write a sorting algorithm
	print "Median g-r:", median(gmr)
	#c)
	i = data.field('i')
	rmi = r-i
	p.plot(rmi,gmr,'.',color='black')
	p.axis([-.5,2.0,-.5,2.0])
	p.xlabel('r-i')
	p.ylabel('g-r')
	p.title('g-r vs. r-i')
	p.savefig('ri_gr_plot')
	p.close()
	#d)
	p.hist2d(rmi,gmr,bins=35,range=((-.5,2.0),(-.5,2.0)))
	p.xlabel('r-i')
	p.ylabel('g-r')
	p.title('g-r vs. r-i')
	p.colorbar()
	p.savefig('ri_gr_hist')
	p.close()
	#
	#2.
	#
	hdulist = fits.open('distribution.fits')
	data = hdulist[0].data
	hdr = hdulist[0].header
	hdulist.close()
	x = data[0]
	f = data[1]
	#normalize the pdf
	bin_size = x[2]-x[1]
	area = n.sum(bin_size*f)
	#rescale f
	f = f/area
	#make x finer resolution
	#spline the data
	spl = InterpolatedUnivariateSpline(x, f)
	#decide on the number of new xs
	xs = n.linspace(-30., 30., 500)
	#b)
	#plot the new pdf
	p.plot(xs, spl(xs), 'g.')
	p.xlabel('x')
	p.ylabel('number of objects')
	p.title('#2 Histogram')
	p.savefig('p2_hist')
	p.close()
	#calculate the expectation values
	Ex = integ.simps(spl(xs)*xs,xs)
	Ex2 = integ.simps(spl(xs)*xs**2,xs)
	stddev = n.sqrt(Ex2 - Ex**2)
	print "Mean x:", Ex
	print "Median x:", median(xs)
	print "Standard Deviation of x:", stddev
	#c)
	#The PDF function that randomly selects from number from
	#the pdf is defined above
	d_i = n.zeros(5)
	for i in range(0,5):
		d_i[i] = PDF(spl,xs)
	print "Mean of d_i:", d_i.mean()
	print "Standard Deviation of d_i:", d_i.std()
	p.hist(d_i,bins=60,range=(-30.,30.))
	p.xlabel('x')
	p.ylabel('number of objects')
	p.title('#2c) Histogram')
	p.savefig('2_histogram')
	p.close()
	#d)Generate 100 samples of each with 5 draws from the pdf
	k = n.zeros((100,5))
	for i in range(0,100):
		for j in range(0,5):
			k[i][j] = PDF(spl,xs)
	k_means = n.zeros(100)
	k_std = n.zeros(100)
	for i in range(0,100):
		k_means[i] = n.mean(k[i][:])
		k_std[i] = n.std(k[i][:])
	p.hist(k_means)
	print "Mean of sample means:", n.mean(k_means)
	#e)Try using N=50,500,5000,50000
	d_51, m_51, std_51 = MC(100,50,spl,xs)
	d_52, m_52, std_52 = MC(100,500,spl,xs)
	d_53, m_53, std_53 = MC(100,5000,spl,xs)
	d_54, m_54, std_54 = MC(100,50000,spl,xs)
	n.save('d_51', d_51)
	n.save('d_52', d_52)
	n.save('d_53',d_53)
	n.save('d_54', d_54)
	#
	#3.
	#a)
	#stats.binom()
	n_values = [20, 30, 50]
	b_values = [0.2, 0.6, 0.5]
	linestyles = ['-', '--', ':']
	x = n.arange(-1, 200)
	# plot the distributions
	for i in range(0,3):    
	    # create a binomial distribution
	    dist = stats.binom(n_values[i], b_values[i])
	    p.plot(x, dist.pmf(x), ls=linestyles[i], color='black', 
	    	label='b=%.1f, n=%i' % (b_values[i], n_values[i]), linestyle='steps-mid')

	p.xlim(-0.5, 35)
	p.ylim(0, 0.25)
	p.xlabel('x')
	p.ylabel('p(x)')
	p.title('Binomial PDF')
	p.legend()
	p.show()
	p.savefig('binomial_pdf')
	p.close()
	#stats.norm
	sigma_values = [0.5, 1.0, 2.0]
	linestyles = ['-', '--', ':']
	mu_values = [-2.0,2.0,0.0]
	x = n.linspace(-10, 10, 1000)

	for i in range(0,3):
	    # create a gaussian / normal distribution
	    dist = stats.norm(mu_values[i], sigma_values[i])
	    p.plot(x, dist.pdf(x), ls=linestyles[i], color='black', label='$\mu=%.1f, \sigma=%.1f$' % (mu_values[i], sigma_values[i]))

	p.xlim(-5, 5)
	p.ylim(0, 0.85)
	p.xlabel('x')
	p.ylabel('p(x)')
	p.title('Gaussian PDF')
	p.legend()
	p.show()
	p.savefig('normal_pdf')
	p.close()
	#stats.poisson
	mu_values = [5, 10, 20]
	linestyles = ['-', '--', ':']
	for i in range(0,3):
		dist = stats.poisson(mu_values[i])
		x = n.arange(-1, 200)
		p.plot(x, dist.pmf(x), ls=linestyles[i], color='black', label='$\mu=%i$' % mu_values[i], linestyle='steps-mid')

	p.xlim(-0.5, 30)
	p.ylim(0, 0.4)
	p.xlabel('x')
	p.ylabel('p(x)')
	p.title('Poisson PDF')
	p.legend()
	p.show()
	p.savefig('poisson_pdf')
	p.close()
	#b)CDFs
	#Binomial CDF
	n_values = [20, 30, 50]
	b_values = [0.2, 0.6, 0.5]
	linestyles = ['-', '--', ':']
	x = n.arange(-1, 200)
	# plot the distributions
	for i in range(0,3):    
	    # create a binomial distribution
	    dist = stats.binom(n_values[i], b_values[i])
	    p.plot(x, dist.cdf(x), ls=linestyles[i], color='black', 
	    	label='b=%.1f, n=%i' % (b_values[i], n_values[i]), linestyle='steps-mid')

	p.xlim(-0.5, 35)
	p.ylim(0, 1.05)
	p.xlabel('x')
	p.ylabel('P(x)')
	p.title('Binomial CDF')
	p.legend(loc=4)
	p.show()
	p.savefig('binomial_cdf')
	p.close()
	#Gaussian CDF
	sigma_values = [0.5, 1.0, 2.0]
	linestyles = ['-', '--', ':']
	mu_values = [-2.0,2.0,0.0]
	x = n.linspace(-10, 10, 1000)

	for i in range(0,3):
	    # create a gaussian / normal distribution
	    dist = stats.norm(mu_values[i], sigma_values[i])
	    p.plot(x, dist.cdf(x), ls=linestyles[i], color='black', label='$\mu=%.1f, \sigma=%.1f$' % (mu_values[i], sigma_values[i]))

	p.xlim(-6, 6)
	p.ylim(0, 1.05)
	p.xlabel('x')
	p.ylabel('P(x)')
	p.title('Gaussian CDF')
	p.legend(loc=4)
	p.show()
	p.savefig('normal_cdf')
	p.close()
	#Poisson CDF
	mu_values = [5, 10, 20]
	linestyles = ['-', '--', ':']
	for i in range(0,3):
		dist = stats.poisson(mu_values[i])
		x = n.arange(-1, 200)
		p.plot(x, dist.cdf(x), ls=linestyles[i], color='black', label='$\mu=%i$' % mu_values[i], linestyle='steps-mid')

	p.xlim(-0.5, 30)
	p.ylim(0, 1.05)
	p.xlabel('x')
	p.ylabel('P(x)')
	p.title('Poisson CDF')
	p.legend(loc=4)
	p.show()
	p.savefig('poisson_cdf')
	p.close()
if __name__ == '__main__':
    main()