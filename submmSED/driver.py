#!/usr/bin/env python
# encoding: utf-8
"""
driver.py

Created by Andrew H. Jaffe on 2010-05-22.
Copyright (c) 2010 Imperial College London. All rights reserved.
"""

from __future__ import division

import sys
import os

import operator

import numpy as np
import matplotlib.pyplot as plt

import MCMC
import likelihood
import data
import model
import getdist


#### pattern after test_binnedCl.py and/or BeamFit/driver.py BeamFit/MAXIPOLBeamData.py

fname_DLC = "./submmSED.txt"
fname_ERCSC = "./submmSED/ercsc_iifscz.txt"

def main(filename=fname_ERCSC, i=0, rotateParams=False, onecomponent=True, getNorm=True, start=None, sigmas=None, 
         nMC=(10000,10000), nDerived=None, noPlots=False, DLC=False, fig0=0, savefig=False, retMCMC=True):
        
        
    ret = []
        
    if not operator.isSequenceType(i):
        idata = [i]
    else:
        idata = i
        
    ## read the data
    if DLC:
        alldata = data.readfluxes_DLC(filename)
    else:
        alldata = data.readfluxes_ERCSC_TopCat(filename)
        
    for i in idata:
        
        try:
            dat = alldata[i]
            name = dat.name + "(z=%4.2f)" % (dat.z)
        except TypeError:
            dat = [alldata[j] for j in i]
            name = " + ".join(str(int(d.name)) for d in dat)
        
        print "Object[s] %s" % name
    
        ## initialize the model (with a class, not an object)
        if getNorm:
            if onecomponent:
                mod = model.submmModel1_normalized_logA
            else:
                mod = model.submmModel2_normalized_logA
            like = likelihood.SEDLikelihood_normalized(data=dat, model=mod)        
        else:
            if onecomponent:
                mod = model.submmModel1
                like = likelihood.SEDLikelihood1(data=dat, model=mod)
            else:
                mod = model.submmModel2            
                like = likelihood.SEDLikelihood2(data=dat, model=mod)
        
        if start is None: 
            start_params = np.asarray(mod.startfrom(random=False))
        else:
            start_params = np.asarray(start)
        if sigmas is None:
            prop_sigmas = start_params/4.0
        else:
            prop_sigmas = sigmas
        
        if nDerived is not None:
            like.nDerived = nDerived
    
        mcmc, ana = MCMC.sampler(like, nMC, prop_sigmas, start_params, plotter=None,
                            fac=None, noCorrelations=True, doBlock=True, rotateParams=rotateParams)

        if not noPlots:

            fig = plt.figure(fig0)
            lnLike = []
            maxlnLike, maxLikeParams = getdist.histgrid(mcmc[-1])
            plt.suptitle(name)
    
            if savefig:
                try:
                    fname = name + savefig
                except TypeError:
                    fname = name
                    
                fig.savefig(fname+"_0.png")
                plt.close(fig)
    
            fig=plt.figure(fig0+1)
            params = ana[0]
            meanmod = mod(*params)
            MLmod = mod(*maxLikeParams)
            try:
                meanmod.plot(dat, wavelength=True, logplot=True)
                MLmod.plot(dat, wavelength=True, logplot=True)
                plt.suptitle(name)    
            except AttributeError:
                pass
            if savefig:
                fig.savefig(fname+"_1.png")
                plt.close(fig)
                
            fig0 += 2
    
        if retMCMC:
            ret.append((mcmc, ana, (maxlnLike, maxLikeParams))) 
        else:
            ret.append((ana, (maxlnLike, maxLikeParams))) 
    
    if len(ret) == 1:
        ret = ret[0]
        
    return ret

def many():
    
    idata =[i*50 for i in range(6,29)]
    
    nMC = (15000,100000)

    print "Two-Component beta = 2"
    ret1 = main("submmSED/ercsc_iifscz.txt", getNorm=True, i = idata, 
                start=(1,2.,10,0.1,2.,20), sigmas=(1,0,2, 1, 0, 2), retMCMC=False,
                nMC=nMC, onecomponent=False, fig0=0, savefig="_2comp_b2")

    print "One-Component"
    ret2 = main("submmSED/ercsc_iifscz.txt", getNorm=True, i = idata, 
                start=(1,2.,10), sigmas=(1,2,2), retMCMC=False,
                nMC=nMC, onecomponent=True, fig0=100, savefig="1comp")
                
    print "One-Component beta = 2"
    ret3 = main("submmSED/ercsc_iifscz.txt", getNorm=True, i = idata, 
                start=(1,2.,10), sigmas=(1,0,2), retMCMC=False,
                nMC=nMC, onecomponent=True, fig0=200, savefig="1comp_b2")
                
    return ret1, ret2, ret3
    

def plotter(sampler):
    """
    make a plot from the results of the sampler
    (could probably make this more generic if plotMod is a member fn of the model)
    for now, only plot the *first* dataset (since for this case that's the 'good' data)
    
    TODO: calculate stats for the amplitude here?
    
    NOT USED IN ABOVE YET
    """

    mod = sampler.like.model
    data = sampler.like.data

    params = where(sampler.prop.sigmas>0)[0]


    ntot = len(sampler.samples)
    stride = 1 #max(1,ntot//sampler.naccept)
    ana = MCMC.chain_analyze(sampler.samples[(ntot//5)::stride,:], params=params)
    vals = sampler.like.model.package(ana[0])
    sigs = sampler.like.model.package(ana[1])

    ### need to take shape into account
    vals = mod.bandpowers(vals)
    sigs = mod.bandpowers(sigs)


    #vals = vals*mod.ellctr*(mod.ellctr+1)/(2*math.pi)
    #sigs = sigs*mod.ellctr*(mod.ellctr+1)/(2*math.pi)


    print vals  #sampler.like.model.fmtstring % tuple(ana[0])
    print sigs  #sampler.like.model.fmtstring % tuple(ana[1])
    print ana[2]


    #m = mod(vals); c1 = m()[0]; ell = N.arange(c1.size)
    #c = c1*ell*(ell+1)/(2*N.pi)
    #pylab.plot(ell, c)


if __name__ == '__main__':
    main()

