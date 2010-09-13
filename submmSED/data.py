from __future__ import division

from numpy import asarray, float64, log, where, nonzero, math, loadtxt, genfromtxt, concatenate
import numpy as np
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

    def __init__(self, freq, flux, sigma, name, z, nu_obs=None, z_alt=None, name_alt=None):
        GaussianData.__init__(self, flux, sigma)        
        self.freq=asarray(freq)
        #self.name = asarray(name)
        self.name = str(name)
        self.z = asarray(z)
        self.freq_obs = nu_obs
        self.z_alt = z_alt
        self.name_alt = name_alt


    def plot(self, lab=True, wavelength=False, fmt=None, logplot=True):
        """ plot the data """
        if wavelength:
            x = speed_of_light/self.freq
            xlab = "wavelength [microns]"
        else:
            x = self.freq
            xlab = "frequency [GHz]"
                        
        if logplot:
#            plt.semilogy(x, self.d, fmt)
            plt.loglog(x, self.d, fmt)
        else:
            plt.plot(x, self.d, fmt)
        one_e = 1.0 - 1e-10
        sigs = [np.min([one_e*self.d, self.sig],axis=0), self.sig]
        plt.errorbar(x, self.d, sigs, fmt=None)
        if lab:
            plt.xlabel(xlab)
            plt.ylabel("flux")


def readfluxes_DLC(filename):
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
        name = str(int(name))
        data.append(submmData(nu_rest, flux, sig, name, z, nu_obs=nu_obs))
        
    return data
        
        
def readfluxes_ERCSC_TopCat(filename, upperlim=2.0, delete_upperlim=False):
    """read fluxes from a TopCat file (complicated format)
    if upperlim, then convert detections at less than (upperlim*sigma) 
    into upper limits
       (specifically, set <f>=0,
            flux<0 gets sig<=2*sig, and flux>0 get sig<=2*flux
    if delete_upperlim, then delete detections at less than (upperlim*sigma) 
    """
    
    # dtype = [('Flux217', np.float64), ('Err217', np.float64), 
    #           ('Flux353', np.float64), ('Err353', np.float64), 
    #           ('Flux545', np.float64), ('Err545', np.float64), 
    #           ('Flux857', np.float64), ('Err857', np.float64), 
    #           ('IRASName', 'S12'),
    #           ('RA', np.float64), ('DEC', np.float64), 
    #           ('Flux12u', np.float64), 
    #           ('Flux25u', np.float64), 
    #           ('Flux60u', np.float64), 
    #           ('Flux100u', np.float64),
    #           ('Zspec', np.float64), ('Z', np.float64),
    #           ('NEDName', 'S12'),
    #           ('Separation', np.float64)]
    #  
    
    errfrac = 0.1   ### fractional error for IRAS
    lambda_IRAS = asarray([12.0, 25.0, 60.0, 100.0]) ## microns
    err_IRAS = asarray([0, 1.0, 0.1, 0.1])
    nu_Planck = asarray([217., 353., 545., 857.])  ## GHz
    
    Planck_idx = (0,2,4,6)
    IRAS_idx = (11,12,13,14)
    
    IRAS_idx = IRAS_idx[1:]   ### ignore 11=12 micron!
    err_IRAS = err_IRAS[1:]
    lambda_IRAS = lambda_IRAS[1:]
    
    nu_IRAS = speed_of_light/lambda_IRAS  ## GHz
    nu_obs = concatenate((nu_Planck, nu_IRAS))


    with open(filename) as f:
        while True:
            l = f.readline()   ### this construction is necessary to avoid 
                               ### "ValueError: Mixing iteration and read methods would lose data "
            if l[0:2] == "+-": break
            
        l = f.readline()  ## delimiter line
        l = f.readline()
        names = [s.strip() for s in l.split("|")]
        
        lines = genfromtxt(f, delimiter="|", comments="+-", dtype=None, usecols=range(1,20))
        
    data = []
    for obj in lines:
        
        name = obj[8].strip()
        name_alt = obj[17].strip()
        zspec = obj[15]
        z = obj[16]

        ## convert Planck fluxes to mJy from Jy
        flux = asarray([1e-3*obj[i] for i in Planck_idx] + [obj[i] for i in IRAS_idx])
        sig  = asarray([1e-3*obj[i+1] for i in Planck_idx] + [ef*obj[i] for i, ef in zip(IRAS_idx, err_IRAS)])

        if delete_upperlim:
            idx_good = np.logical_and(flux>0, flux/sig>upperlim)
            flux = flux[idx_good]
            sig = sig[idx_good]
            nu_obs = nu_obs[idx_good]
            if not np.all(idx_good):
                name += "D"
        elif upperlim:
            idx_lt0 = flux<0
            idx_gt0 = np.logical_and(flux > 0, flux/sig<upperlim)
            idx = np.logical_or(idx_lt0, idx_gt0)
            if np.any(idx):
                flux_gt0 = flux[idx_gt0]
                sig[idx_lt0] = 2*sig[idx_lt0]
                flux[idx] = 0.01*sig[idx]
                sig[idx_gt0] = 2*flux_gt0
                name += "U"

        nu_rest = nu_obs*(1.+zspec)

        data.append(submmData(nu_rest, flux, sig, name, zspec, 
                              nu_obs=nu_obs, z_alt=z, name_alt=name_alt))

    return data
