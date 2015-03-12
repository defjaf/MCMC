from __future__ import division

##### PROBLEM

####### THE FOLLOWING IS NO LONGER RELAVENT
#### self.mod is the *class* of the model, but self.mod.proposal is an actual instance
#### this means that there is only one instance of the proposal per python run!!!
###### could be fixed by re-initing self.proposal in the model constructor, but this is much slower
###### better solution would be a factory function or metaclass which creates a *class* with the appropriate package/unpackage

###### currently, self.mod.proposal is indeed a class, but the package and unpackage methods are set by a factory function.

### 10 June 2010: added derived parameters -- likelihood needs getDerived(*params), nDerived, derivedTexNames
###   should just keep them in the parameter samples array or treat separately????

### Dec 2010: Allow randomized starting points. 
##            currently, setStart just returns if given a bad starting point, which fails later on.
##            should allow some sort of loop?

### AHJ 09/13 -- if random, randomize the restart position (by a certain number of sigma?)
###              see newMCMC
###  TODO      multiple chains -- gelman-rubin convergence test?


import math

### this is actually a forked branch from https://github.com/fnoble/python-progressbar/tree/ipython_support
# from progressbar_ipython import ProgressBar
### `pip install progressbar-ipython` but shows up as just `import progressbar`
from progressbar import ProgressBar

from numpy import (array, float64, exp, log, concatenate, zeros, bool8,
                   reshape, sqrt, identity, where, asarray, ones, int32,
                   arange)
from numpy.random import uniform, seed
from numpy.linalg import LinAlgError
from numpy import empty, isnan, isfinite, isneginf, isinf
import numpy as np

from Likelihood import ZeroPosterior

### assumes likelihood > 0 (so use lnLike) but prior can be 0 (so call
### directly from model)

#debugf = file("debug.out", "w")

from contextlib import contextmanager

@contextmanager
def finishing(obj):
    try:
        yield obj
    finally:
        obj.finish()

class MCMC(object):
    """    Metropolis sampler    """
    
    def __init__(self, like, proposal=None, startProposal=None,
                 nMC=0, startParams=None, doBlock=False):
        """
        initialize the sampler with particular likelihood/posterior class
        i.e., like= Likelihood.Likelihood(data, model)
            must have method lnLike(parameters)
                      member model, with static methods
                                 prior(), package(), unpackage()
        
        if proposal is None, then like.model.proposal must exist
        initialize the proposal with values from startProposal (e.g., sigmas)
        can do further initialization with methods self.prop
        
        if nMC>0 & starting parameters are supplied, start the MCMC
        
        allow 'fixed' parameters via sigma=0 proposal density
        (affects calculation of normalized covariance matrix)
        
        allow updating parameters in blocks:
            doBlock = True
            model class must define nBlock, paramBlocks
                    = [int for each param giving block (0...nBlock-1)]
            save acceptance stats per block!
                        
            TODO:   problem with single params with sig=0 (or NAN)?
                    problem with blocks with all params with sig==0.
                      (latter fixed by explicitly leaving these parameters out of
                      the blocks -- could be automated by checking for sig=0?)
        """
        
        self.like = like
        
        self.doBlock = doBlock
        self.iBlock = 0
     
     # "promote" some of the attributes of these members to members of this class
        
        self.prob = like.lnLike # i.e. the __call__ method in this case...
                                #      or set to lnPost?
        self.mod = self.like.model
     
     #  this way these get pickled with the rest of the instance
        try:
            self.texNames = self.mod.texNames
            self.paramBlocks = self.mod.paramBlocks
            self.nBlock = self.mod.nBlock
            self.fmtstring = self.mod.fmtstring
        except AttributeError:
            pass
        
        self.data = like.data
        
        seed(seed=None)
        
        if proposal is None:
            self.prop = self.mod.proposal()  ### nb. proposal is a class so this initializes it.
        else:
            self.prop = proposal()
        
        if startProposal is not None:
            self.prop.setSigmas(startProposal)
        
        if startParams is not None:
            try:
                self.setStart(startParams)
            except ZeroPriorError:
                print "Can't start at", startParams, ": prior = 0"
                return
            except ZeroPosterior:
                print "Can't start at", startParams, ": posterior = 0"
                return
            
            if nMC>0:
                self.MC_append(nMC)
    
    
    ### TODO: set random from caller!
    def setStart(self, params, random=True, verbose=True):
        """ set the starting point for the chain. 
         AHJ 2014/08: randomize the starting point 
        
        """
        paramarr = self.like.model.unpackage(params)
        self.nparams = len(paramarr)
        self.samples = empty(dtype=float64, shape=(1, self.nparams))
        self.lnPr = empty(dtype=float64, shape=(1,))
        self.accepted = empty(dtype=float64, shape=(1, ))
        prior = self.like.model.prior(*params)
        
        ### AHJ 08/2014: randomize the start by drawing from the proposal
        if random:
            prior = -1 ## force randomization
            i_restart = 0
            while prior <= 0:
                newparams = self.prop.getNewParams(params)      ### want to resample all parameters each time, even if doBlock=True
                                                                ### nb. always start at the given params
                paramarr = self.like.model.unpackage(newparams)
                prior = self.like.model.prior(*newparams)
                i_restart += 1
                
            if i_restart>1 and verbose:
                print "Number of prior resamples: %d" % i_restart
                
        params = newparams
        self.samples[0] = paramarr
        self.prev = params
        
        try:
            self.nDerived = self.like.nDerived
        except AttributeError:
            self.nDerived = self.like.nDerived = 0    
        
        
        if prior<=0:  ### should never get here if random=True
            raise ZeroPriorError(params)
        
        self.lnPr[0] = self.prev_lnPr = self.prob(params) + log(prior)
        
        self.naccept = 1
        if self.doBlock:
            self.accept = ones(dtype=int32, shape=self.like.model.nBlock)
        
        self.accepted[0] = True
        
        if self.nDerived:
            self.derived = empty(dtype=float64, shape=(1, self.nDerived))
            self.derived[0] = self.like.getDerived(*params)
    
    
    def getStart(self):
        " return the starting point for this chain "
        return self.like.model.package(self.samples[0])
    
    
    ## could add a progress bar here?
    def MC_append(self, nMC=1):
        """ add nMC more MCMC iterations to the chain. """
        ### store alpha, probability, accepted or not?
        
        ## store the 'unpackaged' samples in the chain for speed
        
        newsamples = empty(shape=(nMC, self.nparams), dtype=float64)
        newlnPr = empty(shape=nMC, dtype=float64)
        if self.nDerived:
            newDerived = empty(shape=(nMC, self.nDerived), dtype=float64)
        newaccepted = empty(shape=nMC, dtype=bool8)

        with finishing(ProgressBar().start()) as pbar:
            for i in xrange(nMC):  ## do with comprehension?
                samples, newlnPr[i], newaccepted[i], derived = self.sample()
                if self.nDerived:
                    newDerived[i] = derived
                newsamples[i] = self.like.model.unpackage(samples)
                pbar.update(100.0*(i+1)/nMC)
        
        self.samples = concatenate((self.samples, newsamples))
        self.lnPr = concatenate((self.lnPr, newlnPr))
        self.accepted = concatenate((self.accepted, newaccepted))
        if self.nDerived:
            self.derived = concatenate((self.derived, newDerived))
    
    
    ## would it be faster (possibly cython-able) to lift this up into MC_append and refactor the whole range(nMC) loop?
    def sample(self):
        """ get a single next sample from the chain; nb in 'native' format """
        
        accepted = False
        if self.doBlock:
            block = where(asarray(self.like.model.paramBlocks) == self.iBlock)[0] # nb. where() returns a tuple
        else:
            block = None
        
        next = self.prop.getNewParams(self.prev, block=block)
        
        try:
            next_lnPr = self.prob(next)
            if isnan(next_lnPr) or isneginf(next_lnPr):    ### AHJ CHECK -- is this really needed or is below better???
                #self.n_NaN += 1
                raise ZeroPosterior
        except ZeroPosterior:
            ## don't accept, no matter what
            prior = 0
        else:
            prior = self.like.model.prior(*next)
            
        
        # working with lnalpha eliminates some over/underflow errors
        if prior > 0:
            next_lnPr += log(prior)
            lnalpha = next_lnPr-self.prev_lnPr
            lnalpha += self.prop.lndensityRatio(self.prev, next)
            
            # if not isfinite(lnalpha):
            #     print 'lnalpha: next_lnPr=%f, prev_lnPr=%f' % (next_lnPr, self.prev_lnPr)
    
    #        print >> debugf, prior, '|', next, '|', next_lnPr,
        
        ### AHJ CHECK: is this quite right? Does it screw up when the *previous* prob was NaN?
        #### OK, since NaN is never saved as prev
        ### lnalpha<0 check short-circuits the exp(lnalpha)?
        #if prior<=0 or (lnalpha < 0 and exp(lnalpha)<uniform(0,1)):
        ## rewrite as not (not check) to catch NaNs if they appear
        if not prior>0 or (not lnalpha >= 0 and not exp(lnalpha)>uniform(0,1)):
            next = self.prev   ### reject: keep the old parameters
            next_lnPr = self.prev_lnPr
        else:
            accepted = True
            self.naccept += 1
            if self.doBlock:
                self.accept[self.iBlock] += 1
    
    #        print >> debugf, accepted
        
        self.prev = next
        self.prev_lnPr = next_lnPr
        
        if self.doBlock:
            self.iBlock = (self.iBlock+1) % self.like.model.nBlock
            
        if self.nDerived:
            derived = self.like.getDerived(*next)
        else:
            derived = None
        
        return next, next_lnPr, accepted, derived
    
    
    ######## duplicates effort; shouldn't recalculate if burn, stride don't change
    def mean(self, burn=0, stride=1):
        """ mean of the chain, with optional burn and stride (thinning)
        """
        
        params = range(self.nparams)
        self.means = [
            self.samples[burn::stride,i].mean() for i in params ]
        return self.like.model.package(self.means)
    
    def stdev(self, burn=0, stride=1):
        """ std deviation of the chain, with optional burn and stride (thinning)
        """
        params = range(self.nparams)
        chain = self.samples[burn::stride]
        #nsamp = chain.shape[0]   ##nsamp = (chain.shape[0]-burn)//stride  WRONG
        self.stdevs = [ chain[:,i].std() for i in params ]
        
        return self.like.model.package(self.stdevs)
    
    
    def covar(self, burn=0, stride=1, unNormalized=False, stdevs=None, params=None):
        """ normalized covariance of the chain, with optional burn and
        stride (thinning)
            only calculate for sequence of params if not None
              (in which case returned matrix is just the appropriate submatrix)
        """
        if params is None: params = range(self.nparams)
        
        ## nb, mean(), stdev() calculated on *all* params, even if 'fixed'
        ##   (so have the same index i as chain[:,i])
        
        means = self.like.model.unpackage(self.mean(burn, stride))
        if stdevs is None:
            stdevs = self.like.model.unpackage(self.stdev(burn, stride))  ###
        chain = self.samples[burn::stride]
        
        self.covars = zeros((self.nparams, self.nparams), dtype=float64)
        
        if not unNormalized:
            for i in params:
                for j in params:
                    self.covars[i,j] = (
                        ((chain[:,i]-means[i])*(chain[:,j]-means[j])).mean()/
                        (stdevs[i]*stdevs[j]))
        else:
            for i in params:
                for j in params:
                    self.covars[i,j] = (
                        ((chain[:,i]-means[i])*(chain[:,j]-means[j])).mean())
                
        return self.covars
    
    
    def newMCMC(self, burn=0, stride=1, fac=None, nMC=None, noCorrelations=False,
                      doBlock=None, rotateParams=False, random=False):
        """
        start a new chain based on the mean and variances of the present one.
        
        starts at the present mean, with stdev and covar given by the
        present samples.  relies on the proposal density taking the
        stdev as a starting 'width', as well as the normalized
        covariance.  The latter is likely only applicable to gaussian
        or gaussian-like densities (and even stdev may not be
        applicable)
        
        if rotateParams, try to rotate to an orthogonal parameter basis, determined by the previous covariance matrix
        [[ TODO: allow setting the rotation explicitly?? ]]
        
        """
        ### AHJ 09/13 -- if random, randomize the restart position (by a certain number of sigma?)
        
        params = where(self.prop.sigmas>0)[0]  # where returns a tuple
        
        newStart = self.mean(burn, stride)
        stdevs = self.stdev(burn, stride)
        
        if random: newStart += np.random.randn(len(stdevs))*stdevs   ### need to check for 0 prior
        
        try:
            covar = self.covar(burn, stride, params=params)
        except ZeroDivisionError:
            print "newMCMC: ZeroDivisionError"
            ### want to allow stdev=0 for fixed parameters!
            pstdev = self.like.model.unpackage(stdevs)
            pstdev[stdevs<=0.0] = 1.e-4
            covar = self.covar(burn, stride, stdevs=pstdev, params=params)
            
        
        npar = self.like.model.nparam
        if fac is None:
            fac = 2.4/sqrt(npar)
        
        if doBlock is None:
            doBlock=self.doBlock
        
        startP = self.like.model.package(fac*array(self.like.model.unpackage(stdevs)))
        newSampler = MCMC(self.like, startProposal=startP, startParams=newStart,
                          doBlock=doBlock)

        if rotateParams or not noCorrelations:
            newSampler.prop.rotateParams = rotateParams
            noCorrelations = False
            try:
                newSampler.prop.setNormalizedMatrix(covar)
            except AttributeError:
            ## only makes sense for gaussian-like proposal densities
                pass
            except LinAlgError:
                print 'uh oh: non-positive matrix'
                newSampler.prop.setNormalizedMatrix(covar+identity(npar)*0.01)
            except:
                print "shouldn't be here; matrix is:"
                print covar
                raise
        
            
        
        if nMC is not None:
            newSampler.MC_append(nMC)
    
    #        print >> debugf, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    #        debugf.flush()
        
        return newSampler
    
    def copyMCMC(self):
        """ Copy self to new: preserve the starting point and the
        proposal, but not the actual chain, if any
        """
        newSampler = MCMC(self.like)
        newSampler.setStart(self.getStart())
        newSampler.prop = self.prop ## nb this is a reference to the original
                                    ##    want a copy?
        return newSampler
    
    def __getstate__(self):
        odict = self.__dict__.copy() # copy the dict since we change it
        ## replace various entries with their names
        odict['prob'] = odict['prob'].__name__
        odict['mod'] = odict['mod'].__name__
        del odict['prop'] # = odict['prop'].__name__
        del odict['like'] #= odict['like'].__name__
        return odict
    
    ## should write a __setstate__ to restore these attributes, somehow!

class ZeroPriorError(Exception):
    
    def __init__(self, *value):
        self.value = value
    
    def __str__(self):
        return repr(self.value)


def chain_analyze(chain, params=None):
    """
    get second-order statistics for an MCMC chain
    
    assumes shape=(number of MCs, number of parameters)
    
    So this should be called on the appropriate slice as in
    chain_analyze(chain[burnin::stride,:])
    
    returns (means, stddev, normalized covariance matrix)
    
    deprecated -- just use MCMC.mean(), MCMC.stdev(), MCMC.covar()
    
    if params is a seq. of ints, only use those in the analysis
    """
    
    n = chain.shape[1]
    
    if params is None: params = range(n)
    np = len(params)
    
    nsamp = chain.shape[0]
    means = [ chain[:,i].mean() for i in range(n) ]
    
    stdevs = [ chain[:,i].std() for i in range(n) ]
    
    covar=zeros( (n,n), dtype=float64 )
    for i in params:
        for j in params:
            covar[i,j] = (((chain[:,i]-means[i])*(chain[:,j]-means[j])).mean()/
                          (stdevs[i]*stdevs[j]))
    
    return means, stdevs, covar

def sampler(like, nMC, prop_sigmas, start_params, plotter=None, fac=None,
            noCorrelations=False, doBlock=False, rotateParams=False, randomrestart=False):
    """
    sample from the likelihood with a series of MCMC runs given by the
    sequence nMC, using the endpoint of one as the start of the
    next. The first instance starts at start_params with the proposal
    initialized with prop_sigma
    
    if a single chain doesn't produce enough samples for
    statistics, append the same amount again.
    
    takes a callable plotter(sampler) for plotting, displaying
    data
    """
    mod = like.model
    data = like.data
    sampler = []
    
    if rotateParams: 
        noCorrelations = False
    noCorr1 = noCorrelations
    
    nresample = 0
    nresample_max = 10
    
    params = where(like.model.unpackage(prop_sigmas)>0)[0]  # where returns a tuple
    
    for isamp, nMC1 in enumerate(nMC):
        if isamp==0:  ## first MCMC has uncorrelated paramters, proposal width given by prop_sigmas 
            new_s = MCMC(like, startProposal=prop_sigmas, nMC=0,
                              startParams=start_params, doBlock=doBlock)
        else:
            new_s = sampler[isamp-1].newMCMC(burn=nMC[isamp-1]//burnfrac,
                                             stride=stride, nMC=0, fac=fac,
                                             noCorrelations=noCorr1, rotateParams=rotateParams,
                                             random=randomrestart)
        
        burnfrac = nMC1
        
        while True:
            new_s.MC_append(nMC1)
            print "done with chain %d. naccept=%d (%f)" % (
                    isamp, new_s.naccept, new_s.naccept/len(new_s.samples))
            if new_s.doBlock:
                print "per block: ", new_s.accept
            
            if new_s.naccept < 10:
                noCorr1 = True
                print "Not using correlations"
            else:
                noCorr1 = noCorrelations
            
            if new_s.naccept == 0:
                continue
            
            ntot = len(new_s.samples)
            
            stride = 1    ### max(1,ntot//new_s.naccept) ### always use stride=1
            try:
                ana = chain_analyze(new_s.samples[(ntot//burnfrac)::stride,:],
                                    params=params)
                means = array(ana[0])[params]
                stdvs = array(ana[1])[params]
            except ZeroDivisionError:
                print 'ZeroDivisionError; resampling'
                continue
            
            #check for very small variance...
            eps = 1.e-9
            with np.errstate(invalid='raise'):
                try:
                    if max(stdvs/means) < eps:
                        print 'small relative variance; resampling'
                        nresample += 1
                        if nresample > nresample_max: 
                            print 'reached max resample limit: breaking out of loop'
                            break
                        continue
                except (ZeroDivisionError, FloatingPointError):
                    idx = abs(means)<1.e-11
                    if max(stdvs[idx]) < eps:
                        print 'small absolute variance; continuing'
                        continue
            
            
            break   # if you get to here, you're done!
        
        sampler.append(new_s)
        
        if plotter is not None: plotter(new_s)
    
    return sampler, ana

