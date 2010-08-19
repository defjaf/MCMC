from __future__ import division

from numpy import asarray, float64, log, where, nonzero, math, loadtxt
from numpy.random import uniform, normal

import matplotlib.pyplot as plt

from GaussianData import GaussianData

# AHJ: modify/extend/subclass for non-uniform noise

##  multiple SEDs -- should just require initialization with a list of data...

speed_of_light = 299792.458 ### micron GHz

class submmData(GaussianData):
    """
    represents SED measurements.
    assume data at given frequencies (freq[]), with fluxes (flux[]) and errors (sigma[])
    with noise variance sigma^2 (should work for array or scalar sigma)
    """

    def __init__(self, freq, flux, sigma, name, z, nu_obs=None):
        GaussianData.__init__(self, flux, sigma)        
        self.freq=asarray(freq)
        self.name = asarray(name)
        self.z = asarray(z)
        if nu_obs is not None:
            self.freq_obs = nu_obs
        else:
            self.freq_obs = None

    def plot(self, lab=True, wavelength=False, fmt=None):
        """ plot the data """
        if wavelength:
            x = speed_of_light/self.freq
            xlab = "wavelength [microns]"
        else:
            x = self.freq
            xlab = "frequency [GHz]"
                        
        plt.semilogy(x, self.d, fmt)
        plt.errorbar(x, self.d, self.sig, fmt=fmt)
        if lab:
            plt.xlabel(xlab)
            plt.ylabel("flux")


def readfluxes(filename):
    """read fluxes from a DLC file: each line is [name f1 e1 f2 e2 f3 e3 f3 e4 z]"""
    lines = loadtxt(filename, skiprows=2)
    
    lambda_obs = asarray([60.0, 100, 450, 850]) ## microns
    nu_obs = speed_of_light/lambda_obs ### GHz
    data = []
    for obj in lines:
        name, f1, e1, f2, e2, f3, e3, f3, e4, z = obj
        nu_rest = nu_obs*(1.0+z)
        flux = obj[1:9:2]
        sig  = obj[2:9:2]
        data.append(submmData(nu_rest, flux, sig, name, z, nu_obs=nu_obs))
        
    return data
        