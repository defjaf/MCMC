


from numpy import asarray, float64, log, where, nonzero, math, loadtxt, genfromtxt, concatenate
import numpy as np
import os
import fnmatch
from numpy.random import uniform, normal
import csv
import re
import pkg_resources

import matplotlib.pyplot as plt

from ..GaussianData import GaussianData

# AHJ: modify/extend/subclass for non-uniform noise

##  multiple SEDs -- should just require initialization with a list of data...

#color corrections
# 857GHz    1.02 
# 545GHz    1.1
# 353GHz    1.13

verbose = False

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
            xlab = "wavelength [$\mu$m]"
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
            plt.xlabel(xlab, size='x-large')
            plt.ylabel("flux density [Jy]", size='x-large')
            
        plt.ylim(bottom=self.d[self.d>0].min()/10.)
        plt.tick_params(labelsize='large')
        
        

def readfluxes_DLC(filename, format=0, delnu=None):
    """read fluxes from a DLC file: each line is 
    
    format 0: [name f1 e1 f2 e2 f3 e3 f3 e4 z]
    format 1: [z f1 e1 f2 e2 f3]
    format 2: [name f1 ... f10 e1 ... e10]   z~0, 


    headers: Fxxx = xxx micron flux
             Fxxx_yyy = xxx micron = yyy GHz flux
    
    """
    
    if format == 0:
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

    elif format == 1:
        lines = loadtxt(filename, skiprows=1)
        
        lambda_obs = asarray([250., 350., 500.])
        nu_obs = speed_of_light/lambda_obs ### GHz
        for i,obj in enumerate(lines):
            z = obj[0]
            nu_rest = nu_obs*(1.0+z)            
            flux = obj[1::2]
            sig  = obj[2::2]
            name = str(i)+"_"+str(z)
            data.append(submmData(nu_rest, flux, sig, name, z, nu_obs=nu_obs))
            
    elif format == 2:
        lines = loadtxt(filename, skiprows=1)

        lambda_obs = asarray([25.0, 60, 100, 250, 350, 350, 500, 550, 850, 1400])
        print("converting 25u to upper limit")
        nu_obs = speed_of_light/lambda_obs ### GHz
        data = []
        nband = len(lambda_obs)
        for obj in lines:
            name = obj[0]
            z = 0
            nu_rest = nu_obs  ## z=0 for these
            flux = obj[1:nband+1]
            sig  = obj[nband+1:]
            sig[0] = flux[0]   ## convert detection to 1-sigma upper limit
            flux[0] = 0.0
            name = str(int(name))
            data.append(submmData(nu_rest, flux, sig, name, z, nu_obs=nu_obs))
            
    elif format == 3:
        data = readfluxes_peel(filename, delnu=delnu)
        

    return data

def toFloat(s):
    return np.float(s) if s!='' else np.nan

def readfluxes_DLC_2014(filename="./dat/herus_phot.csv", UL25=True, getArp220=True,
                        del157=True, is_csv=None):
    """read fluxes from a DLC 2014 CSV file: each line is 
    
    name [obsID] z f1 e1 f2 e2 ... fn en
    headers: Fxxx_Jy = xxx micron flux  (or S_xxx)
             Exxx_Jy = xxx micron error
             
    [obsID] is an optional column -- should just be ignored by all processing
             
    if there is flux but no error, force error = np.nan (to treat as UL)
    
    if UL=True, treat all 25 micron data as upper limit
    
    if del157=True, delete 157mu point iff there is a 160mu point!
    
    if has_ObsID=True, there is an extra column between name and z 
       (otherwise read the header and check for it)
             
    NB!!! the original numbers file had a typo in the Redshift header name! was "Redfshift"!

    """

###   allow any of S_123, F123, F123_Jy
###   allow z for redshift column name
###   err = 0 --> no data

    Npat = r'Name|Source'
    Epat = r'^E\d+_|^E\d+$'   ### regex patterns 
    Fpat = r'^F\d+_|^S\_\d+$|^F\d+$'
    Zpat = r'Redshift|Redfshift|z|Z'

    data = []
    with open(filename) as ofile:
        is_csv = "," in ofile.readline() if is_csv is None else is_csv
        ofile.seek(0)

        if is_csv:
            hreader = csv.reader(ofile)
            headers = np.array(next(hreader))
            arr = hreader
        else:
            headers = np.array(ofile.readline().split())
            assert headers[0][0] == "#"
            headers[0] = headers[0][1:]
            arr = np.loadtxt(ofile, dtype='S')  ### 'S' seems to pick the minimum possible len
            
        ncol = len(headers)
        
        name_column = np.array([i for i,h in enumerate(headers) if re.match(Npat, h)])
        assert len(name_column)==1
        name_column = name_column[0]

        z_column = np.array([i for i,h in enumerate(headers) if re.match(Zpat, h)])
        assert len(z_column)==1
        z_column = z_column[0]

        flux_columns = np.array([i for i,h in enumerate(headers) if re.match(Fpat, h)])
        err_columns = np.array([i for i,h in enumerate(headers) if re.match(Epat, h)])
        ### assert np.all(err_columns == flux_columns+1)### no longer true!
        
        lambda_obs = np.array([float(headers[i].lstrip('FS_').rstrip('_H_Jy')) for i in flux_columns])  ### lambda in um
        lambda_err = np.array([float(headers[i].lstrip('E_').rstrip('_H_Jy')) for i in err_columns])  
        assert(np.all(lambda_obs==lambda_err))   ### check correct order

        for row in arr:  
            
            if row[name_column]=="Arp220":
                row[name_column] = "Arp220-short"
            row = np.asarray(row)
            if np.all(row==''): continue  ## ignore blank lines

            z = np.float(row[z_column])
            dat = list(zip(row[flux_columns], row[err_columns], lambda_obs))
            dat_compressed = np.array([[np.float(f), toFloat(e), l] for f,e,l in dat if f!=''])

            ## remove err=0 entries
            dat_compressed = np.array([[f,e,l] for f,e,l in dat_compressed if e>0])
            
            if del157:
                l157 = np.where(np.round(dat_compressed[:,2])==157.0)[0]
                l160 = np.where(np.round(dat_compressed[:,2])==160.0)[0]
                if len(l157) and len(l160):
                    if verbose: print(row[name_column], "deleting 157 micron")
                    dat_compressed = np.delete(dat_compressed,l157, axis=0)

            flux = np.array(dat_compressed[:,0])
            sig = np.array(dat_compressed[:,1])
            lambda_obs_1 = np.array(dat_compressed[:,2])
            nu_obs = speed_of_light/lambda_obs_1
            if UL25:
                l25 = np.where(np.round(lambda_obs_1)==25)[0]
                if l25:
                    sig[l25] = flux[l25]
                    flux[l25] = 0
                    
            nu_rest = nu_obs*(1+z)

            data.append(submmData(nu_rest, flux, sig, row[name_column], z, nu_obs=nu_obs))
            
    if getArp220:
        z = 0.018
        Arp220_data = pkg_resources.resource_filename("MCMC.submmSED", "dat/Arp220.txt")
        dat = lambda_obs_1, flux, sig = np.loadtxt(Arp220_data, unpack=True, usecols=(0,1,2))
        nu_obs = speed_of_light/lambda_obs_1
        if UL25:
            l25 = np.where(np.round(lambda_obs_1)==25)[0]
            if l25:
                sig[l25] = flux[l25]
                flux[l25] = 0

        nu_rest = nu_obs*(1+z)

        data.append(submmData(nu_rest, flux, sig, "Arp220", z, nu_obs=nu_obs))
        
        z = 0.018

    return data



def readfluxes_MRR(filename, IRAS_ignore=None, Planck_ignore=None, DLC_ul=False, old_MRR=False, next0=True, colorx=[1.02, 1.1, 1.13, 1.0]):
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
         
   NEW VERSION for file ERCSCiifsczbg.dat
       3  glon,glat,next,dist,
                    ^^^^^
       1  3x,2f10.5,8f11.2,2f10.5,i3,f12.3,3x,
                                  ^^^


    Most of the parameters are from the IIFSCz catalogue, for which there is a 
    readme file at http://astro.ic.ac.uk/~mrr/fss/readmefss
    
    ### colorx gives the factor with which to divide the Planck fluxes at 
         857, 545, 353, 217 (usu. only first three matter?)
    
    """
    
    ## field widths (including whitespace separators)
    delims1 = (13, 11, 11, 4) + 4*(10,) + (3,3,3,6) +\
             16*(8,) + (6,14,11,11,4) + \
             4*(11,) + \
             (3,6,11,4,8,10) + 4*(8,) +\
             (8,8,8,11) + (8,7,11) + 14*(7,) + (7,) +\
             2*(10,) + 8*(11,) + 2*(10,)+(15,)+\
             2*(7,) + (24,23,22)
    
    dtype1 = {
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
    
    delims2 = (13, 11, 11, 4) + 4*(10,) + (3,3,3,6) +\
             16*(8,) + (6,14,11,11,4) + \
             4*(11,) + \
             (3,6,11,4,8,10) + 4*(8,) +\
             (8,8,8,11) + (8,7,11) + 14*(7,) + (7,) +\
             2*(10,) + 8*(11,) + 2*(10,)+ (3,) + (15,)+\
             2*(7,) + (24,23,22)
    
    dtype2 = {
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
            'glon', 'glat', 'next', 'dist',
            'nedtp', 'sdsstp', 'nameNED', 'nameSDSS', 'name2MASS'
        ],
        'formats': 
            ['a12'] + ['f']*2 + ['i'] + ['f']*4 + ['i']*4 +
            ['f']*16 + ['i','f','f','f','i'] +
            ['f']*4 +
            ['i','f','f','i','f','f'] + ['f']*4 +
            ['f']*4 + ['f','i','f'] + ['f']*14 +['i'] +
            ['f']*2 +['f']*8 +['f']*2 + ['i'] + ['f'] + 
            ['a6']*2  + ['a23','a22','a22']
    }
    
    
    if old_MRR:
        delims, dtype = delims1, dtype1
    else:
        delims, dtype = delims2, dtype2
    
    errfrac = 0.1   ### fractional error for IRAS high-quality
    if IRAS_ignore is None:
        IRAS_ignore = set()  ### indices 0...3 of IRAS wavelengths to ignore
    if Planck_ignore is None:
        Planck_ignore = set()

    lambda_IRAS = asarray([12.0, 25.0, 60.0, 100.0]) ## microns
    nu_Planck_all = [857., 545., 353., 217.]  ## GHz
    nu_Planck = asarray([nu for i, nu in enumerate(nu_Planck_all) if i not in Planck_ignore])
    if colorx is not None:
        cx = asarray([c for i, c in enumerate(colorx) if i not in Planck_ignore])
    else:
        cx = asarray(1.0)
        
    if np.any(cx != np.asarray([1.0])):
        print("color corrections:", cx)
            
    lines = np.genfromtxt(filename, dtype=dtype, delimiter=delims)
    
    data = []
    for obj in lines:
        if next0 and obj['next']!=0:
            print("Skipping object %s" % obj['nameIRAS'])
            continue ## ignore anything w/o next==0
            
        z = obj['z']
        name = obj['nameIRAS']
        nu_obs = list(nu_Planck)
        flux = [1e-3*obj[i] for i in ['s%d' % int(f) for f in nu_Planck]]
        sig  = [1e-3*obj[i] for i in ['e%d' % int(f) for f in nu_Planck]]
        flux = [fi/cxi for fi,cxi in zip(flux, cx)]
        sig  = [si/cxi for si,cxi in zip(sig, cx)]
        for i,lam in enumerate(lambda_IRAS):
            if i in IRAS_ignore: 
                continue    ## don't include this frequency
                
            nq = obj['nq%d' % (i+1)]
            flx = obj['s%d' % int(lam)]

            if DLC_ul:
                if int(lam)==25:
                    nq=1   ### force upper limit
            
            if flx <= 0:
                print("skipping lamba=%f for object %s" % (lam, name))
                continue
            
            if nq == 1: # upperlimit 
                sg = flx
                flx = 0.0 ##[0.01*flx
            elif nq == 2:# low qual -- sig = 0.5*flux 
                sg = 0.5*flx
            elif nq == 3: # high qual -- sig = 0.1 * flux
                sg = 0.1*flx
            elif nq == 5:  ##    IRAS Large Galaxy Catalogue, equivalent to 3
                sg = 0.1*flx
            else:
                print('got nq=%d at %s' % (nq, name))
                            
            nu_obs += [speed_of_light/lam]
            flux += [flx]
            sig += [sg]
            
        nu_obs = asarray(nu_obs)
        nu_rest = (1+z)*nu_obs
        data.append(submmData(nu_rest, flux, sig, name, z, nu_obs=nu_obs))
        
    print("Using data in nu_obs order [GHz]:", nu_obs)
    print("Using data in lambda_obs order [mu]:", speed_of_light/nu_obs)

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
        
        lines = genfromtxt(f, delimiter="|", comments="+-", dtype=None, usecols=list(range(1,20)))
        
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

def readfluxes_peel(filename,delnu=None):
    ### peel format: # i j 217 Err 353 Err 545 Err 857 Err 1763 Err 1870 Err 3000 Err 4280 Err 5000 Err 12000 Err 12490 Err 25000 Err
    ### assume z=0??
    ### fluxes in Jy (but really want per Sr?)
    ### AHJ Aug 2012: convert to mJy (but not really needed?)
    nu_obs = np.array([217,353,545,857,1763,1870,3000,4280,5000,12000,12490,25000])
    nnu = len(nu_obs)
    
    data = []
    lines = np.loadtxt(filename)
    if delnu is not None:
        didx = np.searchsorted(nu_obs,delnu)
        print("deleting: ", delnu)
        print("at indices: ", didx)
        nu_obs = np.delete(nu_obs, didx)
        
    for i,obj in enumerate(lines):
        z = 0.0
        nu_rest = nu_obs*(1.0+z)            
        flux = obj[2::2]
        sig  = obj[3::2]
        if delnu is not None:
            flux = np.delete(flux, didx)
            sig  = np.delete(sig,  didx)
        name = '_'.join(str(int(c)) for c in [i,obj[0],obj[1]])
        data.append(submmData(nu_rest, flux, sig, name, z, nu_obs=nu_obs))
    return data
            


def readfluxes_M31(filename="M31/M31Flux-v2.dat",delnu=None):
    ### single SED for M31

    nu, flux, err = np.loadtxt(filename, unpack=True)   #GHz, Jy, Jy
    nnu = len(nu)
    z = 0.0
#     flux *= 1e3 ## convert fluxes to mJy from Jy
#     err *= 1e3

    if delnu is not None:
        didx = np.searchsorted(nu,delnu)
        print("deleting: ", delnu)
        print("at indices: ", didx)
        nu = np.delete(nu, didx)
        flux = np.delete(flux, didx)
        err = np.delete(err, didx)

    nu_rest = nu_obs = nu

    return [submmData(nu_rest, flux, err, "M31", z, nu_obs=nu_obs)]


def readfluxes_mortier(dirname,delnu=None):
    """ read files in Angela Mortier's format (directory with individual files) """

    data = []
    for file in fnmatch.filter(os.listdir(dirname), "*.txt"):
        name = file.split('_')[3]
        fname = os.path.join(dirname,file)
        nu, flux, err = np.loadtxt(fname, skiprows=1, unpack=True)
        flux *= 1e3 ## convert fluxes to mJy from Jy
        err *= 1e3
        nu /= 1e9   ## convert frequencies from Hz to GHz
        data.append(submmData(nu, flux, err, name, 0.0))
        
    return data
        
            
        