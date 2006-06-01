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
import numpy

filename = "data_list.txt"


mapdir = 'cmb/misc-data/MAP'
homedir = os.path.expandvars('${HOME}/home')
if not os.path.exists(homedir):
    homedir = os.path.expandvars('${HOME}')
mapdir = '/'.join( (homedir, mapdir) )

def main(nMC=(1000,)):
    #numarray.Error.pushMode(dividebyzero="warn")
    #numarray.Error.pushMode(all="raise")
    numpy.set_printoptions(precision=4, linewidth=150, suppress=True)
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

    if manybins:
#        bins = [ 2, 11, 21, 31, 41, 51, 61, 81, 101, 121, 141, 161, 181, 201,
#                 221, 241, 261, 281, 301, 351, 401, 451, 501, 551, 601,
#                 651, 701, 801, 901, 1001,  len(llCl)-1]
#        bins = [ 2, 21, 41, 61, 101, 141, 181, 221,  261, 301, 401, 501, 601,
#                 701, 801, 1001, len(llCl)-1]
        bins = [ 2, 21, 61, 141, 221, 301, 501, 701, 1001, len(llCl)-1]
        for i, b in enumerate(bins[:-1]):
            bins[i] = (b, bins[i+1])
        bins = bins[:-1]
        npar = len(bins)

        ell = [int((b[0]+b[1])/2) for b in bins]
        ## start at a reasonable model
        if not testshape:
            Clbins = [b for b in ell if b<len(llCl) ]
            start_params = zeros(shape=(npar,), dtype=float64) + 2000.0
            start_params[0:len(Clbins)] = llCl[Clbins]
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

    pylab.figure(0)
    retval = MCMC.sampler(like, nMC, prop_sigmas, start_params, plotter=plotter,
                        fac=None, noCorrelations=True, doBlock=True)
                        
    if testshape:
        pylab.plot([ll[0], ll[-1]], [1, 1], hold='true')
    else:
        pylab.plot(ll, llCl, hold='true')
        
    pylab.figure(1)
    samples = cat([ s.samples for s in retval[0] ])
    for var in samples.transpose(): pylab.plot(var)     
    
    pylab.figure(2)
    getdist.histgrid(retval[0][-1])
    
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
    print vals  #sampler.like.model.fmtstring % tuple(ana[0])
    print sigs  #sampler.like.model.fmtstring % tuple(ana[1])
    print ana[2]
    
    
    ### or replace with mod.plotmod if written...
    pylab.cla()
    pylab.errorbar(mod.ellctr, vals, yerr=sigs)
    

import hotshot, hotshot.stats
def profrun():
    #profile.Profile.bias=1e-5
    #profile.run('test_binnedCl.main(nMC=100)', 'Clprof.out')
    prof = hotshot.Profile("binnedCl.prof")
    profout = prof.runcall(main, nMC=100)
    prof.close()
    return profout
