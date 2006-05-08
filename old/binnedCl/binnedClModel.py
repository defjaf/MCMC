from __future__ import division

import math
from operator import isNumberType, __and__,  isSequenceType
from copy import copy

from numarray import array, arange, Float64, zeros
import Proposal

class binnedClModel(object):
    """
    model for C_l in bins with specified "window functions"
    doesn't really deal with polarization yet.
    """
    
    def __init__(self, Cb):
        """
        set the parameters to the C_bins = Cb
        """
        self.Cb = Cb

    def __call__(self):
        """
        get a full C_l spectrum for the parameters
        currently works for shapefun = scalar or shapefun=l(l+1)Cl shape

        still need to deal with polarization if present
        """
        # can use the window-function methods from the CosmoMC likelihood?
        self.Cl[0] = self.shapefun
        for Cb, rng in zip(self.Cb, self.bins):
            self.Cl[0, rng[0]:rng[1]+1] *= Cb
            
        self.Cl[0] *= self.ellnorm  # convert to C_l from l(l+1)Cl/2pi
        return self.Cl

    ## use (*C_b) so the list gets 'untupled'
    @staticmethod
    def prior(*C_b):
        
        if reduce(__and__, (c > 0 for c in C_b), True):
            return 1.0
        else:
            return 0.0

    @classmethod
    def setBinning(cls, bins, shapefun=None):
        """
        set up the binning and shape functions
        binning is a sequence [(lmin0, lmax0), (lmin1, lmax1), ...]
        
        shapefun should be a sequence D[l] = [l(l+1)C/2pi] for l=0...ell_max
        if shapefun=None, flat in l(l+1)C_l/2pi (with amplitude 1)
        if shapefun=Scalar, flat in l(l+1)C_l/(2pi)=scalar

        polarization?
        """

        cls.bins = bins # should probably allow just a single sequence
                        # of the start of all bins, followed by lmax
        
        if shapefun is None:
            cls.shapefun = 1.0
        else:
            cls.shapefun = shapefun
        cls.lmax = max(rng[1] for rng in bins)
        cls.lmin = min(rng[0] for rng in bins)

        if isSequenceType(shapefun):
            shapefun=shapefun[:cls.lmax+1]
        
        cls.Cl = zeros(shape=(3,cls.lmax+1), type=Float64)

        ell = arange(cls.lmax+1, type=Float64)
        ell[0] = 1.0
        cls.ellnorm = 2*math.pi/(ell*(ell+1.0))
        cls.ellnorm[0] = 0.0

        cls.nparam = len(bins)
        print 'binning: %d... %d' % (cls.lmin, cls.lmax)
        print 'nparam: ', cls.nparam

    ## don't really need package/unpackage for this case since
    ## parameters are pretty naturally a list (but what about polarization?)
    def package(Cb): return Cb
    def unpackage(Cb): return Cb
    
    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

    unpackage=staticmethod(unpackage)  
    package=staticmethod(package)
