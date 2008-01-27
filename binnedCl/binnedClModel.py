from __future__ import division

import math
from operator import isSequenceType

from numpy import asarray, array, arange, float64, zeros, all, empty, isscalar
import Proposal

## add a function to convert from "parameters" to l(l+1)C_l/(2pi)? 
##   this is just __call__???
## or just cls.binnorm[] factors??

## need a function to go from parameters to C_b
## precomputed binnorm might make sense here
## carlo algorithm:
##C_l = q_{b given l} * C_l^shape
#{\cal C}_b = q_b (\sum_l (l+1/2)/l/(l+1)*{\cal C}_l)/(\sum_l (l+1/2)/l/(l+1))
### nb AHJ shape here is defined for D_l = l(l+1)C_l/(2pi)

# need to set fmtstring -- setBinning

### nb -- 

### should a generic model have a function like BP() or bandpowers()
###    which gives "printable" parameters?
### or use package/unpackage?



class binnedClModel(object):
    """
    model for C_l in bins with specified "window functions"
    doesn't really deal with polarization yet.
    
    the parameters are really the q_b factors rather than bandpowers themselves.
    see bandpowers() for conversion from q_b to C_b
    """
    
    def __init__(self, Cb):
        """
        set the parameters to the C_bins = Cb [[bandpowers]] 
        """
        self.Cb = asarray(Cb)

    def __call__(self):
        """
        get a full C_l spectrum for the parameters
        currently works for shapefun = scalar or shapefun=l(l+1)Cl shape
        
        nb. leaves C_l = shape (not 0) for any non-binned values of ell.
        
        still need to deal with polarization if present
        """
        # can use the window-function methods from the CosmoMC likelihood?
        self.Cl[0] = self.shapefun
        for Cb, rng in zip(self.Cb, self.bins):
            ## nb. bins are (beginning, end) so need the +1 for python
            self.Cl[0, rng[0]:rng[1]+1] *= Cb
            
        self.Cl[0] *= self.ellnorm  # convert to C_l from l(l+1)Cl/2pi
        return self.Cl

    @classmethod
    def bandpowers(cls, qb):
        """return the bandpowers corresponding to a set of qb"""
            
        return qb*cls.BPnorm


    def BP(self, qb=None):
        """return the bandpowers corresponding to the current parameters or a set of qb"""
        
        if qb is None:
            qb = self.Cb
    
        return self.bandpowers(qb)
    
    

    ## use (*C_b) so the list gets 'untupled'
    @staticmethod
    def prior(*C_b):
        
        if all(asarray(C_b)>0):  #reduce(__and__, (c > 0 for c in C_b), True):
            return 1
        else:
            return 0

    @classmethod
    def setBinning(cls, bins, shapefun=None, doBlock=True):
        """
        set up the binning and shape functions
        binning is a sequence [(lmin0, lmax0), (lmin1, lmax1), ...]
        nb. inclusive (non-pythonic), so lmax0<lmin1, etc
        
        shapefun should be a sequence D[l] = [l(l+1)C/2pi] for l=0...ell_max
        if shapefun=None, flat in l(l+1)C_l/2pi (with amplitude 1)
        if shapefun=Scalar, flat in l(l+1)C_l/(2pi)=scalar

        polarization?
        """

        cls.bins = bins # should probably allow just a single sequence
                        # of the start of all bins, followed by lmax
        
        cls.lmax = max(rng[1] for rng in bins)
        cls.lmin = min(rng[0] for rng in bins)
        
        #print 'lmin, lmax, len(shape)=', cls.lmin, cls.lmax, len(shapefun)

        if isSequenceType(shapefun) and len(shapefun)>cls.lmax:
            cls.shapefun=shapefun[:cls.lmax+1]
        elif shapefun is None:
            cls.shapefun = 1.0
        else:
            cls.shapefun = shapefun

        cls.Cl = zeros(shape=(3,cls.lmax+1), dtype=float64)

        ell = arange(cls.lmax+1, dtype=float64)
        ell[0] = 1.0
        cls.ellnorm = 2*math.pi/(ell*(ell+1.0))
        cls.ellnorm[0] = 0.0

        ## bin centres
        cls.ellctr = asarray([(b[0]+b[1])/2 for b in bins])

        cls.nparam = len(bins)
        if doBlock:
            cls.paramBlocks = arange(cls.nparam)
            cls.nBlock = cls.nparam
        else:
            cls.paramBlocks = None
            cls.nBlock = 0
            

        ## conversion factors for qb->bandpowers
        ## I_l[l(l+1)C_l[shape]]/I_l[1] where "logarithmic integral" is
        #       I_l[f_l]=\sum_l[f_l(l+1/2)/l/(l+1)]
        # (nb. shapefun is l(l+1)C_l[shape]), 
        #     so flat in l(l+1)Cl gives 1
        if isscalar(shapefun):
            cls.BPnorm = cls.shapefun
        else:
            cls.BPnorm = bin_spectrum(bins, shapefun)
            
        print 'binning: %d... %d' % (cls.lmin, cls.lmax)
        print 'BPnorm:', cls.BPnorm
        print 'nparam: ', cls.nparam

    ## don't really need package/unpackage for this case since
    ## parameters are pretty naturally a list (but what about polarization?)
    def package(Cb): return asarray(Cb)
    def unpackage(Cb): return asarray(Cb)
    
    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

    unpackage=staticmethod(unpackage)  
    package=staticmethod(package)
    
# make this a separate function so it can be called out of the model context
def bin_spectrum(bins, llCl):
    """ bin the spectrum Cl[] into bins given by bins (list of tuples)"""
    binned_llCl = zeros(len(bins), dtype=float64)
    for i, bin in enumerate(bins):
        ells = arange(bin[0], bin[1]+1, dtype=float64)
        llClbin = llCl[bin[0]:bin[1]+1]
        num = (llClbin*(ells+0.5)/ells/(ells+1.0)).sum()
        denom = ((ells+0.5)/ells/(ells+1.0)).sum()
        binned_llCl[i] = num/denom
    return binned_llCl