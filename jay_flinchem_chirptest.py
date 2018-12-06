# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 19:02:41 2018

@author: eli
"""
from scipy.misc import comb
import matplotlib.pyplot as plt
import numpy as np


def jay_flinchem_chirptest(): 
    from scipy import special
    c1 = 3.5
    c2 = 5.5 
    c3 = 0.0002
    c4 = 6.75*2.*np.pi
    omega = np.array([1.,2.,3.,4])
    gamma = np.array([40.,40.,30.,90.])*2.*np.pi/360.
    
    t = np.linspace(-42, 8., 2500)
    tnorm = t*2.*np.pi
    ai, aip, bi, bip = special.airy(-c3*np.square(tnorm-c4))
    Qr = c1 + c2*ai

    
    A =np.array([0.5,1.,0.25,0.1])
    Aj0 = A*1.
    Aj1 = np.array([0.4,0.4,0.4,0.4])*.97
    print("Aj1 {} ".format(Aj1))
    
    D = np.zeros_like(t)
    

    for i in range(4):
        j = i+1
        phij = gamma[i]*np.sqrt(Qr-1)
        #print("phi {}".format(phij))
        print(Aj0[i])
        Aj = Aj0[i]*(1.-Aj1[i]*np.sqrt(Qr))
        print(Aj)
        D += Aj*np.cos(tnorm*omega[i]-phij)
        
    #zeta = c1 + c2 * 
    plt.plot(t, Qr, 'r', label='Qr')
    plt.plot(t, D + Qr, 'b', label='D')    
    plt.ylim(1, 7.0)
    plt.grid()
    plt.legend(loc='upper left')
    plt.show()


if __name__ == "__main__":
    test2()