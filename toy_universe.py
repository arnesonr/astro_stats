import numpy as n
from matplotlib import pyplot as p
import scipy.stats as stats
#Do this for 10^6 galaxies
theta = n.zeros(1e6)
phi = n.zeros(1e6)
r = n.zeros(1e6)
S_u = n.copy(r) #uniform flux
S_pl = n.copy(S_u) #power law flux
for i in range(0,n.size(r)):
	#randomly distribute theta, phi, r
	theta[i] = n.random.uniform()*2.0*n.pi
	phi[i] = n.arccos(n.random.uniform()*2.0-1.0)
	r[i] = n.random.uniform()**(1./3.)
	S_u[i] = 10./r[i]**2.
#bin a histogram using 4pi/3*(r2^3-r1^3)
hist, bin_edges = n.histogram(r,bins=50)
#normalize hist by the volume of each bin
dV = n.zeros(n.size(hist))
for i in range(0,n.size(hist)):
	dV[i] = ((4.*n.pi)/3.)*(bin_edges[i+1]**3.-bin_edges[i]**3.)
p.plot(hist/dV)
#generate power-law distribution of L n/L ~ L^-2
for i in range(0,n.size(r)):
	S_pl[i] = (10.*n.random.uniform())**2.

p.hist(S_pl,bins=100)

a = n.where(S_u > 1/n.max(r)**2)
S1 = S_u[a]

a = n.where(S_pl > 30.)
S2 = S_pl[a]