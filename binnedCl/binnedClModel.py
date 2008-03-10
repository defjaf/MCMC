from __future__ import division

import math
from operator import isSequenceType

import pylab 

from numpy import (asarray, array, arange, float64, zeros, all, empty, 
                  isscalar, max, dot, min, concatenate)
from numpy import linalg as la

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

### maybe set up a mapping index-list from a 1-d sequence to 2-d TT TE EE

positive_corr = True

class binnedClModel(object):
    """
    model for C_l in bins with specified "window functions"
    doesn't really deal with polarization yet.
    
    the parameters are really the q_b factors rather than bandpowers themselves.
    see bandpowers() for conversion from q_b to C_b
    
    nb. a single vector with all of the values
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
        
        nb. leaves C_l = shape (not 0) for any non-binned values of ell.   ???
        
        still need to deal with polarization if present
        """
        
        ##  loop over = i = {TT, TE, EE [, BB]}...  have class or instance-level has_pol for speed?
        
        # can use the window-function methods from the CosmoMC likelihood?
        ibin = 0
        for iCl in range(self.nCl):
            self.Cl[iCl] = self.shapefun[iCl].copy()
            Cb = self.Cb[ibin:ibin+self.nbins[iCl]]
            for Cbin, rng in zip(Cb,self.bins[iCl]):
                ## nb. bins are (beginning, end) so need the +1 for python
                self.Cl[iCl, rng[0]:rng[1]+1] *= Cbin
            ibin += self.nbins[iCl]
            
            self.Cl[iCl] *= self.ellnorm  # convert to C_l from l(l+1)Cl/2pi
            
        return self.Cl

    @classmethod
    def bandpowers(cls, qb):
        """return list of arrays of the bandpowers corresponding to a set of qb
        
           *** assume qb is a FLAT sequence of the parameters *** 
               -- could use some package/unpackage here???
        """
            
        ret = []
        ibin = 0
        for i in range(cls.nCl):
            ret.append(qb[ibin:ibin+cls.nbins[i]]*cls.BPnorm[i])
            ibin += cls.nbins[i]
        return ret
        
        #return [qb[nbins[]]*cls.BPnorm[i] for i in range(cls.nCl)]
        

    def BP(self, qb=None):
        """return the bandpowers corresponding to the current parameters or a set of qb"""
        
        if qb is None:
            qb = self.Cb
    
        return self.bandpowers(qb)
    
    @classmethod
    def ClCovar(cls, covar):
        """
        convert a <qb qb'> covariance to <Cb Cb'> using the shapefun
        """
        norm = asarray(concatenate(cls.BPnorm))
        return asarray(covar)*norm*norm.reshape(1,cls.nparam)
        

    ## use (*C_b) so the list gets 'untupled'
    ## could enforce positivity but then will need to know which are cross-spectra
    if positive_corr:
        @classmethod
        def prior(cls,*C_b):
            Cb = asarray(C_b)
            if any(Cb[cls.Cltype!=1]<0):   ## Cltype==1 is <TE>
                return 0
            else:
                return 1
            
    else:
        @staticmethod
        def prior(*C_b):
        
            return 1
        
        # if all(asarray(C_b)>0):  #reduce(__and__, (c > 0 for c in C_b), True):
        #     return 1
        # else:
        #     return 0

    @classmethod
    def setBinning(cls, bins, shapefun=None, doBlock=True, nCl=None):
        """
        set up the binning and shape functions
        binning is a sequence [(lmin0, lmax0), (lmin1, lmax1), ...]
        nb. inclusive (non-pythonic), so lmax0<lmin1, etc
        
        shapefun should be a sequence D[l] = [l(l+1)C/2pi] for l=0...ell_max
        if shapefun=None, flat in l(l+1)C_l/2pi (with amplitude 1)
        if shapefun=Scalar, flat in l(l+1)C_l/(2pi)=scalar

        polarization: shapefun is a sequence of D[l]: shapefune[I]=D^I[l]
                      I = TT, TE, EE, BB in that order
        """
        
        if nCl is None:
            if len(bins)==len(shapefun):
                nCl = len(bins)
        if nCl>3:
            raise Exception("Error in shapes of bins, shapefun, giving nCl=%d" % nCl)

        #if nCl==1:
        #    bins = [bins]
            # try:
            #     if bins.ndim==2: 
            #         bins = [bins]
            # except AttributeError:
            #     pass
                    
        cls.nCl = nCl
                
        ## so bins is now bins[iCl][ibin, 0:1]
        print 'nCl:', nCl
        print 'bins:'
        print bins

        cls.bins = bins # should probably allow just a single sequence
                        # of the start of all bins, followed by lmax
        
        
        
        ## gotta be a better way to do this...
        lmax = []
        lmin = []
        for bin in bins:
            lmax.append(array(bin).max())
            lmin.append(array(bin).min())
        cls.lmax=max(lmax)
        cls.lmin=min(lmin)
        
        #cls.lmax = max([b[1] for b in x for x in bins])
        #cls.lmin = min([b[0] for b in x for x in bins])
        # cls.lmax = max(rng[1] for rng in b for b in bins)
        # cls.lmin = min(rng[0] for rng in b for b in bins)
        #cls.lmax = max(rng[1] for rng in b for b in bin for bin in bins)
        #cls.lmax = min(rng[0] for rng in b for b in bin for bin in bins)

        
        print 'lmin, lmax, len(shape)=', cls.lmin, cls.lmax, len(shapefun)

        if isSequenceType(shapefun) and len(shapefun)>cls.lmax:
            cls.shapefun=shapefun[:cls.lmax+1]
        elif isSequenceType(shapefun) and len(shapefun[0])>cls.lmax:
            cls.shapefun=shapefun[:,:cls.lmax+1]
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
        cls.ellctr = [array(b).mean(axis=1) for b in bins]
        
        #cls.ellctr = [asarray([(b[0]+b[1])/2 for b in bin]) for bin in bins]

        cls.nbins = [len(b) for b in bins]
        cls.nparam = sum(cls.nbins)
        
        cls.Cltype = []
        for i, nb in enumerate(cls.nbins): 
            cls.Cltype.extend([i]*nb)
        cls.Cltype = array(cls.Cltype)
        
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
            cls.BPnorm = [bin_spectrum(b, s) for (b, s) in zip(bins, shapefun)]
            
        print 'binning: %d... %d' % (cls.lmin, cls.lmax)
        print 'BPnorm:', cls.BPnorm
        print 'nparam: ', cls.nparam

    ## don't really need package/unpackage for this case since
    ## parameters are pretty naturally a list (but what about polarization?)
    ### for polz, this is probably more efficient, but needs to deal with the different structures in __call__()
    def package(Cb): return asarray(Cb)
    def unpackage(Cb): return asarray(Cb)
    
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
    
    
    
    ### should all of these live here or somewhere else?
    
def plotcorrmat(corrmat, bins=None):
    """ plot a correlation matrix, with a colorbar, optionally with lines to mark different sets of bins"""
    pylab.matshow(corrmat)
    ntot = corrmat.shape[0]-0.5
    if bins is not None:
        xy = bins[0]-0.5
        for b in bins[1:]:
            pylab.plot([-0.5,ntot],[xy, xy])
            pylab.plot([xy, xy],[-0.5,ntot])
            xy+=b
        pylab.xlim(-0.5, ntot-0.5)
        pylab.ylim(ntot-0.5, -0.5)
    pylab.colorbar()
    
    
def orthobin(Cb, corrmat):
    """ 
         rebin bandpowers C_b with corr. mat. <dCb dCb'>=M_bb' 
         to D_p s.t. <dDp dDp'> = diag(sig_p)
         
         nb. Cb not qb
             should be full corr'n matrix not normalized matrix with M_bb=1
             
        AHJ: how to deal with polarization? want to be separately local in ell and XY <- (T, E, B)**2
    """
    
    v, R = la.eigh(corrmat)
    newbins =  dot(N.diag(N.sqrt(v)), R.T)   ## very ineffecient?
    newcorr = dot(R, newbins)
    
    ## make weights out of newbins
    
    ## AHJ: NOT FINISHED
    
    
def fitOffsetLognormal(samples):
    """
    fit an offset lognormal [gaussian in z=ln(C+x)] for <z>, x, sig_z
    
    need to allow for case when C>0 enforced
    """
    
    ## AHJ: NOT FINISHED
    
    
    pass
    
    
    ### TODO NEED TO CONVERT FROM q_B FISHER MATRIX TO C_B MATRIX!!!!!
def FisherWindows(F, bins=None, isCovar=False):
    """ 
        calculate the effective bandpower window functions from the inverse covariance matrix F=C^{-1}
        so WB_l/l = \sum_{l' in B} F_ll'/\sum_{{all l},{l' in B}} F_{ll'}
        
        if bins is not present, just return W_B as a function of bin number, otherwise return full W_Bl
        
    """
    
    if isCovar:
        fish = la.inv(F)
    else:
        fish = F
        
    #Wbb = (fish/fish.sum(axis=0)).T   ## nb. transpose to make each *row* normalized
    Wbb = (fish/fish.diagonal()).T   ## nb. transpose to make each *row* normalized
    # Wbb = fish.copy()
    # for row in Wbb: row /= row.sum()
    
    nbin = Wbb.shape[0]
    
    if bins is None:
        return Wbb
    else:
        ### return W_Bl in three arrays at each B: TT, TE, EE
        
        lmax = []
        for bin in bins: lmax.append(array(bin).max())
        lmax=max(lmax)
        
        print "FisherWindows: lmax=", lmax
        
        ##maxell = max([bin[0][1] for bin in (spec for spec in bins)])
        
        WBl = zeros((nbin, 3, lmax+1), dtype=float64)
        for ibin in xrange(nbin):    ### there must be a more numpyish way to do this...
            jbin = 0
            for ispec, spec in enumerate(bins):
                for bin in spec:
                    WBl[ibin, ispec, bin[0]:bin[1]+1]=Wbb[ibin, jbin]/(bin[1]-bin[0]+1)# * 0.5*(bin[0]+bin[1])
                    jbin += 1
                
        return WBl
