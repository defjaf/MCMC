from __future__ import division

import numarray
from numarray.random_array import uniform

import MCMC
import Likelihood
from BeamFit import BeamModel, BeamData
import math

def main(nMC=1000):
    numarray.Error.pushMode(dividebyzero="warn")
    mod = BeamModel.GaussianBeamModel2D
    sim_model = mod( (0,0), (1,2), math.pi/6 )
    data1 = BeamData.BeamSim(sim_model, 1000, sigma=0.1, amplitude=10.0,
                            xrng=(-1,1), yrng=(-1,1))
    data2 = BeamData.BeamSim(sim_model, 1000, sigma=0.1, amplitude=10.0,
                            xrng=(-1,1), yrng=(-1,1))

    data=[data1, data2]
  #  data = data1

    mod.setxyRange((-1,1))    ## class variables: sets the prior for all instances
    
    like = Likelihood.Likelihood(data=data, model=mod)

    npar = 5
    prop_sigmas = ( (0.1, 0.1), (0.05, 0.05), 0.1 )

    start_params = ( tuple(uniform(-1,1,2)), tuple(uniform(0,1,2)),
                     uniform(0,math.pi/2) )

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
    profile.Profile.bias=1e-5
    profile.run('testMCMC.main()', 'MCMCprof.out')
    

