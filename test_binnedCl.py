
from __future__ import division


import math
import string
import os.path

import pylab
import pyfits
import MCMC
import getdist
from binnedCl.binnedClLikelihood import binnedClLikelihood
from binnedCl.binnedClModel import binnedClModel
from ClData import ClData
from pylab import *
#import numarray
#from numarray.random_array import uniform
#from numarray import arange, array, float64, Error, transpose, zeros, ones
from numpy.random import uniform
from numpy import arange, array, float64, transpose, zeros, ones
from numpy import concatenate as cat
import numpy as N

filename = "data_list.txt"


mapdir = 'cmb/misc-data/MAP'
homedir = os.path.expandvars('${HOME}/home')
if not os.path.exists(homedir):
    homedir = os.path.expandvars('${HOME}')
mapdir = '/'.join( (homedir, mapdir) )

def main(nMC=(1000,), gridPlot=True):
    #numarray.Error.pushMode(dividebyzero="warn")
    #numarray.Error.pushMode(all="raise")
    N.set_printoptions(precision=4, linewidth=150, suppress=True)
    mod = binnedClModel

    #mapf = '/'.join( (mapdir, 'models/cmb_04546021.fits') )
    #mapd = pyfits.getdata(mapf)
    #ll = arange(mapd.shape[0])
    #norm = 1e6
    #ClTT = array(norm**2*mapd.field(0))
    #llCl = ClTT*ll*(ll+1)/(2*math.pi)

    #max_ell = len(llCl)-1
    #max_ell = 2100

    tmp = N.fromfile("CarloClModel.dat", sep=" ")
    tmp.shape = -1, 6
    ells = int_(tmp.T[0])  ## nb doesn't usually start at l=0
    
    max_ell = max(ells)
    ll = arange(max_ell+1)
    llClTT = N.zeros(max_ell+1)
    llClEE = N.zeros(max_ell+1)
    llClTE = N.zeros(max_ell+1)
    llClTT[ells] = tmp.T[1]
    llClEE[ells] = tmp.T[2]
    llClTE[ells] = tmp.T[3]

    manybins = True
    onebin = False

    testshape = True

    data = ClData.getClData(filename, no_pol=True)
    
    # print "unsetting beam and calib uncertainty"
    # for d in data:
    #     d.beam_uncertain=False
    #     d.calib_uncertainty = 0.0
    #     #d.has_xfactors=False

    if manybins:
#        bins = [ 2, 11, 21, 31, 41, 51, 61, 81, 101, 121, 141, 161, 181, 201,
#                 221, 241, 261, 281, 301, 351, 401, 451, 501, 551, 601,
#                 651, 701, 801, 901, 1001,  max_ell]
##        bins = [ 2, 21, 41, 61, 101, 141, 181, 221,  261, 301, 401, 501, 601,
#        bins = [ 45, 61, 101, 141, 181, 221,  261, 301, 401, 501, 601,
#                 701, 801, 1001, max_ell]
        bins = [50, 101, 151, 201, 251, 301, 351, 401, 451, 501, 551, 601, 651, 701, 
                751, 801, 851, 901, 971, 1031, 1091, 1151, 1211, 1271, 1331, 1391, 
                1451, 1511, 1571, 1651, 1751, 1851, 1950, len(llClTT)-1]
##        bins = [ 50, 61, 141, 221, 301, 501, 701, 1001, len(llCl)-1]
##        bins = [ 2, 61, 221, 501, 1001, max_ell]
        for i, b in enumerate(bins[:-1]):
            bins[i] = (b, bins[i+1]-1)  ## nb. non-pythonic: beginning and end
        bins = bins[:-1]
        npar = len(bins)

        ell = [int((b[0]+b[1])/2) for b in bins]
        ## start at a reasonable model
        if not testshape:
            Clbins = [b for b in ell if b<len(llClTT) ]
            start_params = zeros(shape=(npar,), dtype=float64) + 2000.0
            start_params[0:len(Clbins)] = llClTT[Clbins]
            shape = 1.0
            prop_sigmas = zeros(npar, float64) + 100.0
        else:
            shape = llClTT
            start_params = ones(shape=(npar,), dtype=float64)
            prop_sigmas = zeros(shape=(npar,), dtype=float64) + 0.1
            

    elif onebin:
        bins = [(2, len(llClTT)//2), (len(llClTT)//2+1,len(llClTT)-1)]
        npar = len(bins)
        shape = llClTT
        start_params = ones(shape=(npar,), dtype=float64)
        prop_sigmas = zeros(shape=(npar,), dtype=float64) + 0.5

    mod.setBinning(bins, shapefun=shape)

    ell = [(b[0]+b[1])/2 for b in bins]
    print 'bins:'
    print ell
    
    like = binnedClLikelihood(data=data, model=mod)

    pylab.figure(0)
    
    fac = 2.4/sqrt(len(ell))   ## equiv to fac=None
    fac*=3
    
    retval = MCMC.sampler(like, nMC, prop_sigmas, start_params, plotter=plotter,
                        fac=fac, noCorrelations=True, doBlock=True)
                        
    # if testshape:
    #     pylab.plot([ll[0], ll[-1]], [1, 1], hold='true')
    # else:

    pylab.plot(ll, llClTT, hold='true')
        
    pylab.figure(1)
    samples = cat([ s.samples for s in retval[0] ])
    for var in samples.transpose(): pylab.plot(var)
    
    s = retval[0][-1]
    
    if gridPlot:
        pylab.figure(2)
        getdist.histgrid(s)
    else:
        getdist.printvals(s)

    
    params=xrange(s.samples.shape[1])    
    for param in params:
        s1 = s.samples.transpose()[param]

        mean = mod.bandpowers(s1.mean())
        stdv = mod.bandpowers(s1.std())
        print '%d %f %f' % (ell[param], mean[param], stdv[param])        
    
    
    #numarray.Error.popMode()
    return retval

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
    

    ### or replace with mod.plotmod if written...
    pylab.cla()
    pylab.errorbar(mod.ellctr, vals, yerr=sigs, fmt='bo')
    #m = mod(vals); c1 = m()[0]; ell = N.arange(c1.size)
    #c = c1*ell*(ell+1)/(2*N.pi)
    #pylab.plot(ell, c)
    

########################################################
def getlike(ibin=1):
    mod = binnedClModel

    mapf = '/'.join( (mapdir, 'models/cmb_04546021.fits') )
    mapd = pyfits.getdata(mapf)

    ll = arange(mapd.shape[0])
    norm = 1e6
    ClTT = array(norm**2*mapd.field(0))
    llCl = ClTT*ll*(ll+1)/(2*math.pi)

    manybins = True
    onebin = False

    testshape=False

    data = ClData.getClData(filename, no_pol=True)
    
    for d in data:
        d.beam_uncertain=False
        d.calib_uncertainty = 0.0
        #d.has_xfactors=False

    if manybins:
#        bins = [ 2, 11, 21, 31, 41, 51, 61, 81, 101, 121, 141, 161, 181, 201,
#                 221, 241, 261, 281, 301, 351, 401, 451, 501, 551, 601,
#                 651, 701, 801, 901, 1001,  len(llCl)-1]
#        bins = [ 2, 21, 41, 61, 101, 141, 181, 221,  261, 301, 401, 501, 601,
#                 701, 801, 1001, len(llCl)-1]
        #bins = [ 2, 21, 61, 141, 221, 301, 501, 701, 1001, len(llCl)-1]
        bins = [50, 101, 151, 201, 251, 301, 351, 401, 451, 501, 551, 601, 651, 701, 
                751, 801, 851, 901, 971, 1031, 1091, 1151, 1211, 1271, 1331, 1391, 
                1451, 1511, 1571, 1651, 1751, 1851, 1950]
        for i, b in enumerate(bins[:-1]):
            bins[i] = (b, bins[i+1]-1)  ## nb. non-pythonic: beginning and end
        bins = bins[:-1]
        npar = len(bins)

        ell = [int((b[0]+b[1])/2) for b in bins]
        ## start at a reasonable model
        if not testshape:
            Clbins = [b for b in ell if b<len(llCl) ]
            start_params = zeros(shape=(npar,), dtype=float64) + 2000.0
            start_params[0:len(Clbins)] = llCl[Clbins]
            start_params[where(start_params==0)]=2000.0
            shape = 1.0
            prop_sigmas = zeros(npar, float64) + 100.0
        else:
            shape = llCl
            start_params = ones(shape=(npar,), dtype=float64)
            prop_sigmas = zeros(shape=(npar,), dtype=float64) + 0.1
            

    elif onebin:
        bins = [(2, len(llCl)//2), (len(llCl)//2+1,len(llCl)-1)]
        npar = len(bins)
        shape = llCl
        start_params = ones(shape=(npar,), dtype=float64)
        prop_sigmas = zeros(shape=(npar,), dtype=float64) + 0.5

    mod.setBinning(bins, shapefun=shape)

    ell = [(b[0]+b[1])/2 for b in bins]
    print 'bins:'
    print ell
    
    like = binnedClLikelihood(data=data, model=mod)
    
    #### plot the likelihood as a function of power in a single bin
    
    bps = N.linspace(0, 100*start_params[ibin], 100)
    likearr = N.empty_like(bps)
    for i, bp in enumerate(bps):
        pars = start_params
        pars[ibin] = bp
        likearr[i] = like.lnLike(pars)
        
    plot(bps, likearr-max(likearr))
    
    return like, likearr
    
    
########################################################

import hotshot, hotshot.stats
def profrun():
    #profile.Profile.bias=1e-5
    #profile.run('test_binnedCl.main(nMC=100)', 'Clprof.out')
    prof = hotshot.Profile("binnedCl.prof")
    profout = prof.runcall(main, nMC=100)
    prof.close()
    return profout
