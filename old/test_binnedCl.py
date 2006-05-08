from __future__ import division

import math
import string
import os.path

from pylab import *
import pyfits
import MCMC
from binnedCl.binnedClLikelihood import binnedClLikelihood
from binnedCl.binnedClModel import binnedClModel
from ClData import ClData
from pylab import *
import numarray
from numarray.random_array import uniform
from numarray import arange, array, Float64, Error, transpose, zeros, ones


filename = "data_list.txt"


mapdir = 'cmb/misc-data/MAP/'
homedir = os.path.expandvars('${HOME}/home')
if not os.path.exists(homedir):
    homedir = os.path.expandvars('${HOME}')
mapdir = '/'.join( (homedir, mapdir) )

def main(nMC=1000):
    #numarray.Error.pushMode(dividebyzero="warn")
    numarray.Error.pushMode(all="raise")
    mod = binnedClModel

    mapf = '/'.join( (mapdir, 'models/cmb_04546021.fits') )
    mapd = pyfits.getdata(mapf)

    ell = arange(mapd.shape[0])
    norm = 1e6
    ClTT = array(norm**2*mapd.field(0))
    llCl = ClTT*ell*(ell+1)/(2*math.pi)

    manybins = True
    onebin = False

    data = ClData.getClData(filename, no_pol=True)

    if manybins:
        bins = [ 2, 11, 21, 31, 41, 51, 61, 81, 101, 121, 141, 161, 181, 201,
                 221, 241, 261, 281, 301, 351, 401, 451, 501, 551, 601,
                 651, 701, 801, 901, 1001, 2001]
        for i, b in enumerate(bins[:-1]):
            bins[i] = (b, bins[i+1])
        bins = bins[:-1]
        npar = len(bins)

        ell = [int((b[0]+b[1])/2) for b in bins]
        Clbins = [b for b in ell if b<len(llCl) ]
        start_params = zeros(npar, Float64) + 2000.0
        ## start at a reasonable model
        start_params[0:len(Clbins)] = llCl[Clbins]
        shape = 1.0
        prop_sigmas = zeros(npar, Float64) + 100.0

    elif onebin:
        bins = [(2, len(llCl)//2), (len(llCl)//2+1,len(llCl)-1)]
        npar = len(bins)
        shape = llCl
        start_params = ones(shape=(npar,), type=Float64)
        prop_sigmas = zeros(npar, Float64) + 0.5

    mod.setBinning(bins, shapefun=shape)

    ell = [(b[0]+b[1])/2 for b in bins]
    print 'bins:'
    print ell
    
    like = binnedClLikelihood(data=data, model=mod)


    sampler = MCMC.MCMC(like, startProposal=prop_sigmas, nMC=nMC,
                        startParams=start_params)

    print "done with first chain. naccept=", sampler.naccept

    stride = nMC//sampler.naccept
    ana = MCMC.chain_analyze(sampler.samples[(nMC//5)::stride,:])
    vals = sampler.like.model.package(ana[0])
    sigs = sampler.like.model.package(ana[1])
    print vals
    print sigs
    print ana[2]

    try:
        sampler2 = sampler.newMCMC(burn=nMC//5, stride=stride, nMC=nMC, fac=0.1)
    except:
        retval = sampler
    else:
        print "done with second chain. naccept=", sampler2.naccept
        
        stride = nMC//sampler2.naccept
        an2 = MCMC.chain_analyze(sampler2.samples[(nMC//10)::stride,:])
        val2 = sampler2.like.model.package(an2[0])
        sig2 = sampler2.like.model.package(an2[1])
    
        print val2
        print sig2
        print an2[2]
        
        retval = (sampler, sampler2)

    numarray.Error.popMode()
    return retval

import hotshot, hotshot.stats
def profrun():
    #profile.Profile.bias=1e-5
    #profile.run('test_binnedCl.main(nMC=100)', 'Clprof.out')
    prof = hotshot.Profile("binnedCl.prof")
    profout = prof.runcall(main, nMC=100)
    prof.close()
    return profout
