from __future__ import division

import math
import string
import os.path

import pyfits
import MCMC
import binnedClLikelihood
import binnedClModel
import ClData
from pylab import *
import numarray
from numarray.random_array import uniform
from numarray import arange, array, Float64, Error, transpose, zeros


filename = "data_list.txt"

def main(nMC=1000):
    numarray.Error.pushMode(dividebyzero="warn")
    mod = binnedClModel.binnedClModel

    data = []
    for dataset in file(filename):
        if dataset[0] != "#":
            print string.strip(dataset)
            set = ClData.ClData(string.strip(dataset))
            print "Got %s" % (set.name)
            data.append(set)

    bins = [ 2, 11, 21, 31, 41, 51, 61, 81, 101, 121, 141, 161, 181, 201, 221,
             241, 261, 281, 301, 351, 401, 451, 501, 551, 601, 651, 701, 801,
             901, 1001]
    for i, b in enumerate(bins[:-2]):
        bins[i] = (bins[i], bins[i+1])
    bins = bins[:-2]
    
    mod.setbinning(bins, shape=1000.0)    ## class variables: sets the prior for all instances
    
    like = binnedClLikelihood.binnedClLikelihood(data=data, model=mod)

    npar = len(bins)
    prop_sigmas = zeros(npar, Float64) + 10.0

    start_params = zeros(npar, Float64) + 1000.0

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
    
    sampler2 = sampler.newMCMC(burn=nMC//5, stride=stride, nMC=nMC, fac=0.1)
    
    print "done with second chain. naccept=", sampler2.naccept

    stride = nMC//sampler2.naccept
    an2 = MCMC.chain_analyze(sampler2.samples[(nMC//10)::stride,:])
    val2 = sampler2.like.model.package(an2[0])
    sig2 = sampler2.like.model.package(an2[1])
    
    print val2
    print sig2
    print an2[2]
        
    numarray.Error.popMode()
    return sampler, sampler2

import profile
def profrun():
    profile.Profile.bias=5e-6
    profile.run('testMCMC.main()', 'MCMCprof.out')
    

