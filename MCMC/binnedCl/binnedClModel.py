from __future__ import division

import math
from operator import isSequenceType

import matplotlib.pyplot as plt

import pylab

from numpy import (asarray, array, arange, float64, zeros, all, empty, 
                  isscalar, max, dot, min, concatenate, log, pi, inf, abs, sqrt)
from numpy import linalg as la

import numpy as N

import scipy.optimize as So
import scipy.special as Ss

from .. import Proposal

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
        return asarray(covar)*norm*norm.reshape(cls.nparam,1)
        

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
            
        print 'nparam: ', cls.nparam
        #print 'BPnorm:', cls.BPnorm

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
    plt.matshow(corrmat)
    ntot = corrmat.shape[0]-0.5
    if bins is not None:
        xy = bins[0]-0.5
        for b in bins[1:]:
            plt.axvline(xy)
            plt.axhline(xy)
            #plt.plot([-0.5,ntot],[xy, xy])
            #plt.plot([xy, xy],[-0.5,ntot])
            xy+=b
        plt.xlim(-0.5, ntot-0.5)
        plt.ylim(ntot-0.5, -0.5)
    plt.colorbar()
    
    
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

    
### HKE: fit to cumulative distributions
### (and more generally don't use particular functional form, but smoothed "best" transformation to cumulative Gaussian) 
### to get cumulative dist:
###    s = samples[i].sort(); plot(s, range(len(s)))
###    h = histogram(samples, bins=10000); plot(h[1], h[0].cumsum())
### so histogram(samples)[0].cumsum() gives the actual (unnormalized) distribution over 10000 equally-spaced points
### or samples.sort() gives the positions for the unnormalized dist given by range(len(samples))
        
    
class oln(object):
    """offset lognormal likelihood
       this is really just to package it up so it has its own namespace
    """
    
    def __init__(self, C):
        self.C = C

    def chi2(self, zbar, sigz2, x):
        return ((log(self.C+x)/zbar-1)**2).mean()*zbar**2/sigz2
    
    @staticmethod
    def like1(par, C):
        (zbar, sigz2, x) = par
        return log(2*pi*sigz2)+(log(C+x)/zbar-1)**2*zbar**2/sigz2 + 2*log(C+x)
    
    def like(self, par):
        """ return -2lnP(par|samples)"""
        (zbar, sigz2, x) = par
        ### need to check that sigz2>0, C+x>0?
        
        #### HAVE I LEFT OFF THE 1/(C+x) factor????
        
        nsamp = len(self.C)
        if sigz2<0 or x+min(self.C)<0: 
            # raise unphys(par)  ## no way to catch this in So.fmin()
            return inf
        else:
            return log(2*pi*sigz2)+self.chi2(zbar, sigz2, x) + 2*(log(self.C+x)).mean()  ## last term from 1/(C+x)
    
    def derivs(self, par):
        ### need to effect of ln(C+x) term
        (zbar, sigz2, x) = par
        dL_dsigz2 = (1.0/sigz2)*(1.0-self.chi2(zbar, sigz2, x))
        dL_dzbar = (2.0/sigz2)*(zbar - (log(self.C+x)).mean())
        dL_dx = (2.0/sigz2)*((log(self.C+x)-zbar)/(self.C+x)).mean() + (2.0/(self.C+x)).mean()  ## last term from 1/(C+x)
        return array([dL_dsigz2, dL_dzbar, dL_dx])

    def cum(self, par, bins=None):
        """ cumululative OLN """
        
        (zbar, sigz2, x) = par
        
        if bins is None:
            bins = N.sort(self.C)
        elif N.isscalar(bins):
            bins = N.arange(self.C.min(), self.C.max(),bins)
        
        # like = log(2*pi*sigz2)+(log(bins+x)/zbar-1)**2*zbar**2/sigz2 + 2*log(bins+x)  ## last term from 1/(C+x)
        # like = exp(-0.5*like)
        
        ### now need cum sum of like * (delta bin)
        
        ## analytic form!
        
        return 0.5*(Ss.erf((zbar-log(x))/sqrt(2*sigz2))-Ss.erf((zbar-log(bins+x))/sqrt(2*sigz2)))
        
    def KSnorm(self, par):
        return max(abs(self.cum(par)-N.linspace(0,1,len(self.C))))


def fitOffsetLognormal_cum(samples, full_output=0, do_plot=1): 
    """minimize the difference between the actual cumulative sample distribution and the offset lognormal"""
    o = oln(samples)

    #x_0 = max(0,-1.1*min(samples))
    x_0 = 1.1*abs(min(samples))
    zbar_0 = (log(samples+x_0)).mean()
    sigz2_0 = ((log(samples+x_0)-zbar_0)**2).mean()
    par_0 = array([zbar_0, sigz2_0, x_0])
        
    print 'Starting values:', par_0
    print 'Starting chi2:', o.chi2(zbar_0, sigz2_0, x_0)
    print 'Starting KSnorm: ', o.KSnorm(par_0)
    
    f = So.fmin(o.KSnorm, par_0, maxfun=100000, maxiter=100000, xtol=0.00001, ftol=0.00001)
    if do_plot==1:
        bins = N.sort(samples)
        plt.plot(bins, o.cum(f))
        plt.plot(bins, N.linspace(0,1,len(samples)))
    elif do_plot==2:
        bins = N.linspace(samples.min(), samples.max(),100)
        o2 = oln(bins)
        plt.plot(bins,N.exp(-0.5*o2.like1(f, bins)))
        hh = N.histogram(samples, bins=bins, normed=True)
        plt.plot(hh[1], hh[0], ls='steps' )
        
    zbar_f, sigz2_f, x_f = f
    print 'Final chi2:', o.chi2(zbar_f, sigz2_f, x_f)
    print 'Final KSnorm: ', o.KSnorm(f)

    return f
    
fitOffsetLognormal = fitOffsetLognormal_cum
    
def fitOffsetLognormal_like(samples, full_output=0):
    """
    fit an offset lognormal [gaussian in z=ln(C+x)] for <z>, x, sig_z
    
    need to allow for case when C>0 enforced
    """

    o = oln(samples)
    

    #x_0 = max(0,-1.1*min(samples))
    x_0 = 1.1*abs(min(samples))
    zbar_0 = (log(samples+x_0)).mean()
    sigz2_0 = ((log(samples+x_0)-zbar_0)**2).mean()
    par_0 = array([zbar_0, sigz2_0, x_0])
        
    print 'Starting values:', par_0
    print 'Starting chi2:', o.chi2(zbar_0, sigz2_0, x_0)
    print 'Starting derivs:', o.derivs(par_0)
    
    
    f = So.fmin_l_bfgs_b(o.like, par_0, fprime=o.derivs)
    #f = So.fmin_cg(o.like, par_0, fprime=o.derivs, full_output = full_output)
    #f = So.fmin_cg(o.like, par_0, fprime=None, full_output = full_output)
    #f = So.fmin(o.like, par_0, full_output = full_output, maxfun=100000, maxiter=100000)
    
    return f
        
    
    ### nb. NEED TO CONVERT FROM q_B FISHER MATRIX TO C_B MATRIX!!!!!
def FisherWindows(F, bins=None, isCovar=False):
    """ 
        calculate the effective bandpower window functions from the inverse covariance matrix F=C^{-1}
        so WB_l/l = \sum_{l' in B} F_ll'/\sum_{{all l},{l' in B}} F_{ll'}
        
        if bins is not present, just return W_B as a function of bin number, otherwise return full W_Bl
        
        wbl=0.d0
        do i=1,num_bins
           do j=1,num_bins
              sumfish = fisherbb(i,i)
              do l=lbin(j,1),lbin(j,2)
                 wbl(i,type(j),l) =fisherbb(j,i)/(lbin(j,2)-lbin(j,1))/sumfish
              enddo
           enddo
        enddo
    """
    
    if isCovar:
        fish = la.inv(F)
    else:
        fish = F
        
    nbin = fish.shape[0]
    
    if bins is None:
        
        Wbb = (fish.T/fish.diagonal()).T  ## transpose to apply the same factor to each row.
        return Wbb
        
    else:
        
        ### return W_Bl in three arrays at each B: TT, TE, EE
        
        lmax = []
        for bin in bins: lmax.append(array(bin).max())
        lmax=max(lmax)
                
        ##lmax = max([bin[0][1] for bin in (spec for spec in bins)])
        
        WBl = zeros((nbin, 3, lmax+1), dtype=float64)
        for ibin in xrange(nbin):    ### there must be a more numpyish way to do this...
            jbin = 0
            for ispec, spec in enumerate(bins):
                for bin in spec:
                    #WBl[ibin, ispec, bin[0]:bin[1]+1]=Wbb[ibin, jbin]/(bin[1]-bin[0])
                    WBl[ibin, ispec, bin[0]:bin[1]+1]=fish[ibin, jbin]/(bin[1]-bin[0])/fish[ibin, ibin]
                    jbin += 1
                
        return WBl



class unphys(Exception):
    pass