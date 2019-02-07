# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 19:02:41 2018

@author: eli
"""
from scipy.misc import comb
import matplotlib.pyplot as plt
import numpy as np


def lowpass_gain(lambdaf,m,q = 1.):
    out = q*np.power(1./(2.-2.*np.cos(lambdaf)),float(m))
    return out

def spectral_density(lambdaf,lambdac,rho,nn,q=None):
    n = float(nn)
    rho2 = rho*rho
    num = np.zeros_like(lambdaf)
    for j in range(nn+1):
        for k in range(nn+1):
            jkdiff = float(j-k)
            jksum = float(j+k)
            term = np.power(-1.,jksum)*comb(nn,j)*comb(nn,k)*np.cos(lambdac*jkdiff)*np.cos(lambdaf*jkdiff)*np.power(rho,jksum)
            num += term 
            print("Term {} {} {}".format(j,k,term))
    denom= 1.+ 4.*rho2*np.cos(lambdac)**2. + rho2*rho2 - 4.*rho*(1.+rho2)*np.cos(lambdac)*np.cos(lambdaf)+2.*rho2*np.cos(2.*lambdaf)
    denom = np.power(denom,n)
    if not q is None:
        return q*num/denom
        
    leftfrac = (1.-rho2)**(2.*nn-1.)
    outden = 0.
    for i in range(nn):
        rho_to_i = rho**i
        outden += np.square(comb(nn-1,i)*rho_to_i)

    outden *= 2.*np.pi
    print(outden)
    leftfrac /= outden
    print(leftfrac)
    print(num)
    print(denom)
    return leftfrac*num/denom

# These are the order of the trend (nlow) and cyclical (nn) components
nlow = 2
nn = 4
rho = .975
lambdac = 2.*np.pi/99.

#Frequencies in response plot
nsample = 10000
lambdaf = np.linspace(0.,2.*np.pi/16,nsample)  


# Get the gains of the components of the model, evaluated at the frequencies in lambdaf
# Lambdac represents one cycle per day
sdense0 = lowpass_gain(lambdaf,nlow,q=1.)     # subtidal response
sdense1 = spectral_density(lambdaf,1.*lambdac,rho,nn,q=1.)   #diurnal
sdense2 = spectral_density(lambdaf,2.*lambdac,rho,nn,q=1)    # 2* lambdac => semidiurnal
sdense3 = spectral_density(lambdaf,3.*lambdac,rho,nn,q=1)
sdense4 = spectral_density(lambdaf,4.*lambdac,rho,nn,q=1)
sdense5 = spectral_density(lambdaf,.35*lambdac,rho,nn=nn,q=1) # this is an "extra" subtidal

# Need to update total density if you change above
totaldense = sdense0+sdense1+sdense2+sdense3 +sdense4 +sdense5
sdense0 = sdense0/(totaldense+1.)
sdense1 = sdense1/(totaldense+1.)
sdense2 = sdense2/(totaldense+1.)
sdense3 = sdense3/(totaldense+1.)
sdense4 = sdense4/(totaldense+1.)
sdense5 = sdense5/(totaldense+1.)


lambdaflab = np.linspace(0,6,nsample)
plt.plot(lambdaflab,sdense0,label="sub")
plt.plot(lambdaflab,sdense1,label="D1",linewidth=3)
plt.plot(lambdaflab,sdense2,label="D2")
plt.plot(lambdaflab,sdense3,label="D3")
plt.plot(lambdaflab,sdense4,label="D4")
plt.plot(lambdaflab,sdense5,label="D5")
plt.grid()
plt.legend()
plt.show()

   