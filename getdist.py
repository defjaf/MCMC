"""
analyze the output of an MCMC run; functions generally accept either an
instance of class MCMC.MCMC or a shape=(nsamples, nparams) array of samples
(and possibly a 1-d array of lnPr

"""
from __future__ import division

import sys
import os
import os.path
import copy
import numpy as np
import matplotlib.pyplot as plt

from itertools import islice, groupby


def GaussPDF(x, mean=0, stdev=1):
    return np.exp(-0.5*(x-mean)**2/stdev**2)/np.sqrt(2*np.pi)/stdev

def getallsamples(MCMC, derived=True):
    try:
        s = MCMC.samples
        if derived and MCMC.nDerived:
            s = np.hstack([s,MCMC.derived])
    except AttributeError:
        s = MCMC
    return s

def printvals(MCMC, params=None, lnLike=None, derived=True):
    """ like loop over hist(), but don't plot """
    
    s = getallsamples(MCMC, derived=derived)
        
    try:
        if lnLike is None: lnLike = MCMC.lnPr
    except AttributeError:
        lnLike = None

    if lnLike is not None: 
        maxlnLike = max(lnLike)
        maxLikeParams = s[lnLike.argmax()]
        print 'Max ln likelihood %f at parameters:' % maxlnLike
        print maxLikeParams

    if params is None:
        params=xrange(s.shape[1])
    
    for param,ML in zip(params, maxLikeParams):
        s1 = s.T[param]
        
        mean = s1.mean()
        stdv = s1.std()  ## stddev(s1)
        print 'mean = %f +- %f | ML = %f' % (mean, stdv, ML)        
        
    return (maxlnLike, maxLikeParams)
    

def hist(MCMC, param, nbins=10, gauss=True, orientation='vertical', axis=None, derived=True):
    """
    for the MCMC output (or just the list of samples),
    use pylab to make a histograms of parameter number=param
    if gauss==True, also plot a gaussian with the data's mean and stdev
         (only use the 
    """
    
    s = getallsamples(MCMC, derived=derived)        
        
    s1 = s.transpose()[param]
        
    if axis is None: 
        axis=plt.gca()

    hist = axis.hist(s1,  bins=nbins, orientation=orientation)

    if gauss:
        stdv = s1.std()  ## stddev(s1)
        mean = s1.mean()
        print 'mean = %f +- %f' % (mean, stdv)
        if stdv > 0:
            smin, smax = min(hist[1]), max(hist[1])
            norm = s1.size * (smax-smin)/nbins  # sum(hist[0])
            ss = np.arange(smin, smax, (smax-smin)/5/nbins)

            axis.hold(True)
            axis.plot(ss, norm*GaussPDF(ss, mean, stdv), 'k')
        
    return hist
    
    
    
def hists(MCMC, nbins=30, params=None, orientation='vertical',
          nrow=None, ncol=None, derived=True):
    """ plot histograms of samples """
    
    s = getallsamples(MCMC, derived=derived)        

    if params is None:
        npar = s.shape[1]
        params=xrange(npar)
    else:
        npar=len(params)

    if nrow is None and ncol is None:
        sqrtnpar = sqrt(npar)
        nrow = int(sqrtnpar)
        if sqrtnpar == nrow:
            ncol = nrow
        else:
            ncol = int(npar/nrow+1)
        
    for ip in xrange(npar):
        plt.subplot(nrow, ncol, ip+1)   ### 1...n instead of 0...n-1!
        hist(s, ip, nbins=nbins, orientation=orientation)
    

def scatter2d(MCMC, params, **kwargs):
    """
    make a 2d scatterplot of params[0], params[1] from the MCMC.

    if kwarg lnLike is not None, color-code the points by the value

    TODO: overlay gaussian contours
          color code
          
    """

    derived = kwargs.setdefault('derived', False)
    s = getallsamples(MCMC, derived=derived)
    del kwargs['derived']


    s1 = s.transpose()[params[0]]
    s2 = s.transpose()[params[1]]

    lnLike = kwargs.get('lnLike', None)
    if lnLike is not None:
        lnLike=copy.copy(lnLike)  ### necessary? following modifies lnLike
        lnLike -= max(lnLike)  ## normalize to P=1 at the maximum
        del kwargs['lnLike']
        kwargs['c']=lnLike
    
    if 'axis' in kwargs:
        axis = kwargs['axis']
        del kwargs['axis']
    else:
        axis = plt.gca()

    return axis.scatter(s1,s2, s=.01, edgecolors='none', **kwargs)

def histgrid(MCMC, params=None, nbins=30, labels=None, lnLike=None, quiet=False, derived=True, maxlnLike=[], noPlot=False):
    """
    make a 2d grid of histograms and scatterplots of the selected parameters
    """

    s = getallsamples(MCMC, derived=derived)
        
    try:
        if lnLike is None: lnLike = MCMC.lnPr
    except AttributeError:
        lnLike = None
        
    if lnLike is not None: 
        maxlnLike = max(lnLike)
        maxLikeParams = s[lnLike.argmax()]
        print 'Max ln likelihood %f at parameters:' % maxlnLike
        print maxLikeParams
        
    if noPlot and lnLike is not None:  ### hmmmm, should refactor?
        return (maxlnLike, maxLikeParams)
        
    if params is None:
        # try:   ### removed June 2010 to deal with derived parmaters. Don't need this case?
        #     params = MCMC.paramBlocks
        #     npar = len(params)
        # except AttributeError:            
            npar = s.shape[1]
            params=xrange(npar)
    else:
        npar=len(params)

    if labels is None:
        try:
            labels=MCMC.texNames
            if derived: 
                labels += MCMC.like.derivedTexNames
        except AttributeError:
            pass

    nrow = ncol = npar

    #norm = plt.normalize(vmin=-10, vmax=0)
    norm = plt.normalize()

    fig = plt.gcf()

    for ipar1, par1 in enumerate(islice(params, 1, npar)):        ### rows 
        for ipar2, par2 in enumerate(islice(params, ipar1+1)):          ### columns
                  
            ax=fig.add_subplot(nrow, ncol, npar*ipar1+ipar2+1)
            ax.hold(False)
            #splot = \ 
            scatter2d(MCMC, (par2, par1), lnLike=lnLike, norm=norm, axis=ax, derived=derived)
                    ### par2, par1 since x-axis along columns
            ax.set_xticklabels([])
            if ipar2 != 0:
                ax.set_yticklabels([])
            elif labels is not None:
                ax.set_ylabel(labels[par1])
                ax.set_yticks(ax.get_ylim())
                
                
    for ipar, par in enumerate(params):
        ax=fig.add_subplot(nrow, ncol, npar*(npar-1)+ipar+1)
        ax.hold(False)
        hist(MCMC, par, nbins=nbins, axis=ax)
        if ipar != 0:
            ax.set_yticklabels([])
        else:
            ax.set_yticks(ax.get_ylim())

        if labels is not None:
            ax.set_xlabel(labels[par])
            ax.set_xticks(ax.get_xlim())
            
    if not quiet: plt.draw()
    
    if lnLike is not None:
        return (maxlnLike, maxLikeParams)
    

### no longer needed -- array.std() now uses 1/N
def stddev(arr):
    return arr.std()*np.sqrt((arr.size-1)/arr.size)


def convertSampleFile(filename):
    """ 
    convert a CosmoMC-style output file to arrays for use with the 
    functions in this module
   
    input format (per line):
    
    number_of_repetitions -lnlike [parameter_value[i] for i in nparams]
    """
    fil = open(filename, 'r')

    lnLike = []; samples = []
    for lin in fil:
        cols = lin.split()
        nsamp = int(cols[0])
        npar = len(cols) - 2
        lnLike.extend([-float(cols[1])]*nsamp)
        samples.extend([float(c) for c in cols[2:]]*nsamp)

    fil.close()
    samples = np.array(samples)
    nrow = len(samples)/npar
    samples.shape=(nrow, npar)

    return np.array(lnLike), samples

def main(filename, burn=None, labels=None):

    lnLike, samples = convertSampleFile(filename)

   # if labels is None:
   #     labels = ['EE1', 'EE2', 'EE3', 'BB1', 'BB2', 'BB3']

    #print 'lnLike: ', lnLike.shape, lnLike.dtype
    #print 'samples: ', samples.shape, samples.dtype
    histgrid(samples[burn:], labels = labels)
    plt.show()
   
   
def key(filename):
    return filename.rsplit('_', 1)[0]
   
def doall(dir=None, burn=None, labels=None, thin=None, quiet=False, keyfn=None):
    """
    XX1 is ell=2 to ell=150
    XX2 is ell=151 to ell=693
    XX3 is ell=694 to ell=1999
    
    To get physical C_ell from the numbers in the chain you need to multiply by C_shape(mid-range ell).
    where
    C_shape(ell) = 5089.38 / (ell(ell+1))
    OR           = 1000.0 / (2*ell+1)
    OR            = 1.0
    
    for the different chains.    
    
    For ell=422 (the middle of the non-junk bin) these numbers are:
    2.865e-02,  1.185e+00,   and 1.0 respectively.
    
    TO DO: combine chains from same params
           require C>0
    """

    if dir is None: dir = "/Users/jaffe/Desktop/Downloads/MAXIPOL chains/"
    if labels is None: 
        labels = ["EE1","EE2","EE3","BB1","BB2","BB3","EB1","EB2","EB3"]

    if burn is None: burn=0
    if thin is None: thin=1

    for root, dirs, files in os.walk(dir):
        for key, names in groupby(files, keyfn):
            print 'key:', key
            i=0
            for name in names:
                base, ext = os.path.splitext(os.path.basename(name))
                if ext == '.txt' and len(base)>0:
                    print 'file:', base
                    lnLike1, samples1 = convertSampleFile(dir+name)
                    if i==0:
                        lnLike=lnLike1[burn::thin]; samples=samples1[burn::thin]
                    else:
                        lnLike=np.concatenate((lnLike, lnLike1[burn::thin]))
                        samples=np.concatenate((samples, samples1[burn::thin]))
                    i+=1
                    
            if i>0:
                plt.clf()
                histgrid(samples, labels=labels, lnLike=-lnLike, quiet=quiet)
                plt.savefig(dir+key+'.png')
                
               
if __name__=="__main__":
    main(sys.argv[1])
