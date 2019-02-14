#!/usr/bin/env python
# encoding: utf-8
"""
reprocess_greybody.py

Created by Andrew H. Jaffe on 2011-03-30.
Copyright (c) 2011 Imperial College London. All rights reserved.

driver.py in the current directory does not correctly display the grey-body fluxes in the last columns of the output file
SHOULD GO BACK AND FIX THIS

"""

import sys
import os

import numpy as np

from .model import totalflux

speedOfLight = 299792.  ## micron GHz
nu2, nu1 = speedOfLight/8.0, speedOfLight/1000.0 ### rest-frame microns -> GHz


# calculate the fluxes and flux ratios. can work for different models 
#but would have to change the dtype and the various tab['...'] uses
def rp(fname, npar=4, nt=2):
    
    
    ### requires fixing the lines that start with " F" and those with "A(z="
    dt = np.dtype([
        ('name', 'S12'),
        ('zz', 'S10'),
        ('z', np.float),
        ('ML', np.float, (npar,)), 
        ('mean', np.float, (npar,)), 
        ('sig', np.float, (npar,)), 
        ('dlnLike', np.float), 
        ('evidence', np.float, (2,)), 
        ('dat', np.float, (ndat,2)),
        ('flux', np.float, (nt,))
    ])

    tab = np.loadtxt(fname, skiprows=1, dtype=dt)
    
    beta = 2.0
    logA_mean = tab['mean'][:,(0,2)]
    T_mean = np.array(tab['mean'][:,(1,3)], dtype=np.float)
    logA_ML = tab['ML'][:,(0,2)]
    T_ML = np.array(tab['ML'][:,(1,3)], dtype=np.float)
    
    flux_ML = 10.0**logA_ML * totalflux(beta, T_ML, nu1, nu2)
    flux_mean = 10.0**logA_mean * totalflux(beta, T_mean, nu1, nu2)
    
    ratio_ML = flux_ML[:,0]/flux_ML[:,1]
    ratio_mean = flux_mean[:,0]/flux_mean[:,1]
    
    return tab, flux_ML, flux_mean, ratio_ML, ratio_mean


coldnames = [
    "F11211+1805",
    "F02285-4315",
    "F01380-2909",
    "F14366+0534",
    "F08564-0442",
    "F03505-3730",
    "F22260-5840"
    ]
    
def getratios(coldnames=coldnames):
    tab, flux_ML, flux_mean, ratio_ML, ratio_mean = rp("figs_Planckfinal/next0_model0.npy")
    
    print(('name', 'mean logA0', 'mean T0', 'mean logA1', 'mean T1','ML logA0', 'ML T0', 'ML logA1', 'ML T1', 
           'mean flux0', 'mean flux1',  'ML flux0', 'ML flux1', 'mean ratio', 'ML ratio'))
    for c in coldnames:
        idx = np.where(tab['name']==c)[0]
        lin = np.hstack( ((c,), tab['mean'][idx][0], tab['ML'][idx][0],flux_mean[idx][0], flux_ML[idx][0], ratio_mean[idx][0], ratio_ML[idx][0]) )
        print(lin)
    

if __name__ == '__main__':
    pass