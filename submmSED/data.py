from __future__ import division

from numpy import asarray, float64, log, where, nonzero, math, loadtxt
from numpy.random import uniform, normal

import matplotlib.pyplot as plt

from GaussianData import GaussianData

# AHJ: modify/extend/subclass for non-uniform noise

##  multiple SEDs -- should just require initialization with a list of data...

class submmData(GaussianData):
    """
    represents SED measurements.
    assume data at given frequencies (freq[]), with fluxes (flux[]) and errors (sigma[])
    with noise variance sigma^2 (should work for array or scalar sigma)
    """

    def __init__(self, freq, flux, sigma, name, z):
        GaussianData.__init__(self, flux, sigma)        
        self.freq=asarray(freq)
        self.name = asarray(name)
        self.z = asarray(z)

    def plot(self, lab=True):
        """ plot the data """
        plt.semilogy(self.freq, self.d)
        plt.errorbar(self.freq, self.d, self.sig)
        if lab:
            plt.xlabel("frequency")
            plt.ylabel("flux")


def readfluxes(filename):
    """read fluxes from a DLC file: each line is [name f1 e1 f2 e2 f3 e3 f3 e4 z]"""
    lines = loadtxt(filename, skiprows=2)
    
    nu = asarray([60.0, 100, 450, 850]) ## GHz
    data = []
    for obj in lines:
        name, f1, e1, f2, e2, f3, e3, f3, e4, z = obj
        flux = obj[1:9:2]
        sig  = obj[2:9:2]
        data.append(submmData(nu, flux, sig, obj[0], obj[-1]))
        
    return data
        