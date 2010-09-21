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
        GaussianData.__init__(self, asarray(flux), asarray(sigma))        
        self.freq=asarray(freq)
        #self.name = asarray(name)
        self.name = str(name)
        self.z = asarray(z)
        self.freq_obs = asarray(nu_obs)
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
        #sigs = [np.min([one_e*self.d, self.sig],axis=0), self.sig]
        sigs = self.sig
        plt.errorbar(x, self.d, sigs, fmt=fmt)
        if lab:
            plt.xlabel(xlab)
            plt.ylabel("flux")
            
        plt.ylim(ymin=self.d[self.d>0].min()/10.)


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




def readfluxes_MRR(filename):
    """
    read fluxes from an MRR file. 
    
           read(4,135)nameIRAS, ra, dec, posFlag,
         1 s12, s25, s60, s100, nq1, nq2, nq3, nq4,
         2 am1, am2, am3, am4, am5, am6, am7, am8, em1, em2, em3, em4,
         3 em5,em6, em7,em8, photFlag,FINT, EFINT,
         4 zspec, zspecFlag, zneur, zneurerr, ztem,z,
         5 j2, av1, err1, zneurFlag, amb2, alb,
         6 alp1, alp2, alp3, alp4,alcirr, alsb, ala220, alagn, alir,
         7 nirtem, errir3,als12, als25, als60, als90, als100, als110,
         8 als140,als160, als250, als350,als500, als850,als1250,
         1  als1380,nirflag,
         2  ra1,dec1,s857,e857,s217,e217,s353,e353,s545,e545,
         3  glon,glat,dist,
         9 nedtp, sdsstp, nameNED, nameSDSS, name2MASS

    135    format(a12,1x,2(f10.5,1x),i3,1x,4(f9.3,1x),4(i2,1x),3x,
         1 16(f7.2,1x),i3,3x,f13.5,1x,f10.5,1x,f10.6,1x,i3,1x,
         2 4(f10.6,1x),
         2 i2,1x,f5.2,1x,f10.3,1x,i3,1x,f7.2,1x,f7.2,3x,4(f7.4,1x),
         3 4(f7.2,1x),3x,f7.2,1x,i4,3x,f10.3,1x,14(f6.2,1x),i3,1x,
         1  3x,2f10.5,8f11.2,2f10.5,f12.3,3x,
         4 2(a6,1x),a23,1x,a22,1x,a22)


    Most of the parameters are from the IIFSCz catalogue, for which there is a 
    readme file at http://astro.ic.ac.uk/~mrr/fss/readmefss
    
    """
    
    ## field widths (including whitespace separators)
    delims = (13, 11, 11, 4) + 4*(10,) + (3,3,3,6) +\
             16*(8,) + (6,14,11,11,4) + \
             4*(11,) + \
             (3,6,11,4,8,10) + 4*(8,) +\
             (8,8,8,11) + (8,7,11) + 14*(7,) + (7,) +\
             2*(10,) + 8*(11,) + 2*(10,)+(15,)+\
             2*(7,) + (24,23,22)
    
    dtype = {
        'names': [
            'nameIRAS', 'ra', 'dec', 'posFlag',
            's12', 's25', 's60', 's100', 'nq1', 'nq2', 'nq3', 'nq4',
            'am1', 'am2', 'am3', 'am4', 'am5', 'am6', 'am7', 'am8', 'em1', 'em2', 'em3', 'em4',
            'em5', 'em6', 'em7', 'em8', 'photFlag', 'FINT', 'EFINT',
            'zspec', 'zspecFlag', 'zneur', 'zneurerr', 'ztem', 'z',
            'j2', 'av1', 'err1', 'zneurFlag', 'amb2', 'alb',
            'alp1', 'alp2', 'alp3', 'alp4', 'alcirr', 'alsb', 'ala220', 'alagn', 'alir',
            'nirtem', 'errir3', 'als12', 'als25', 'als60', 'als90', 'als100', 'als110',
            'als140', 'als160', 'als250', 'als350', 'als500', 'als850', 'als1250',
            'als1380', 'nirflag',
            'ra1', 'dec1', 's857', 'e857', 's217', 'e217', 's353', 'e353', 's545', 'e545',
            'glon', 'glat', 'dist',
            'nedtp', 'sdsstp', 'nameNED', 'nameSDSS', 'name2MASS'
        ],
        'formats': 
            ['a12'] + ['f']*2 + ['i'] + ['f']*4 + ['i']*4 +
            ['f']*16 + ['i','f','f','f','i'] +
            ['f']*4 +
            ['i','f','f','i','f','f'] + ['f']*4 +
            ['f']*4 + ['f','i','f'] + ['f']*14 +['i'] +
            ['f']*2 +['f']*8 +['f']*2 + ['f'] + 
            ['a6']*2  + ['a23','a22','a22']
    }
    
    errfrac = 0.1   ### fractional error for IRAS high-quality
    IRAS_ignore=set()  ### indices 0...3 of IRAS wavelengths to ignore

    lambda_IRAS = asarray([12.0, 25.0, 60.0, 100.0]) ## microns
    nu_Planck = asarray([857., 545., 353., 217.])  ## GHz
        
    lines = np.genfromtxt(filename, dtype=dtype, delimiter=delims)
    
    data = []
    for obj in lines:
        z = obj['z']
        name = obj['nameIRAS']
        nu_obs = list(nu_Planck)
        flux =  [1e-3*obj[i] for i in ['s%d' % int(f) for f in nu_Planck]] 
        sig  =  [1e-3*obj[i] for i in ['e%d' % int(f) for f in nu_Planck]] 
        for i,lam in enumerate(lambda_IRAS):
            if i in IRAS_ignore: 
                continue
                
            nu_obs += [speed_of_light/lam]
            nq = obj['nq%d' % (i+1)]
            flx = obj['s%d' % int(lam)]
            if nq == 1: # upperlimit 
                flux += [0.0] ##[0.01*flx]
                sig += [flx]
            elif nq == 2:# low qual -- sig = 0.5*flux 
                flux += [flx]
                sig += [0.5*flx]
            elif nq == 3: # high qual -- sig = 0.1 * flux
                flux += [flx]
                sig += [0.1*flx]
            elif nq == 5:  ##    IRAS Large Galaxy Catalogue, equivalent to 3
                flux += [flx]
                sig += [0.1*flx]
            else:
                print 'got nq=%d at %s' % (nq, name)
                
        nu_obs = asarray(nu_obs)
        nu_rest = (1+z)*nu_obs
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
                flux[idx] = 0 ##0.01*sig[idx]
                sig[idx_gt0] = 2*flux_gt0
                name += "U"

        nu_rest = nu_obs*(1.+zspec)

        data.append(submmData(nu_rest, flux, sig, name, zspec, 
                              nu_obs=nu_obs, z_alt=z, name_alt=name_alt))

    return data
