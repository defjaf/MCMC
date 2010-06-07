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

import numpy as np

import MCMC

import likelihood
import data
import model

import getdist


#### pattern after test_binnedCl.py and/or BeamFit/driver.py BeamFit/MAXIPOLBeamData.py


def main(filename="./submmSED.txt", i=0, rotateParams=False):
    
    
    nMC = 10000,10000
    
    
    ## read the data
    alldata = data.readfluxes(filename)
    dat = alldata[i]
    
    ## initialize the model (with a class, not an object)
    mod = model.submmModel2

    start_params = np.asarray(mod.startfrom(random=False))
    prop_sigmas = start_params/4.0

    like = likelihood.SEDLikelihood2(data=dat, model=mod)
    
    mcmc, ana = MCMC.sampler(like, nMC, prop_sigmas, start_params, plotter=None,
                        fac=None, noCorrelations=True, doBlock=True, rotateParams=rotateParams)

    getdist.histgrid(mcmc)
    
    return mcmc, ana


def plotter(sampler):
    """
    make a plot from the results of the sampler
    (could probably make this more generic if plotMod is a member fn of the model)
    for now, only plot the *first* dataset (since for this case that's the 'good' data)
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

