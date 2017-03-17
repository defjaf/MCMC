"""
analyze the output of an MCMC run in file form
"""


import sys
import numpy
import pylab
from itertools import islice


def GaussPDF(x, mean=0, stdev=1):
    return numpy.exp(-0.5*(x-mean)**2/stdev**2)/numpy.sqrt(2*numpy.pi)/stdev


def hist(MCMC, param, nbins=10, gauss=True, orientation='vertical'):
    """
    for the MCMC output (or just the list of samples),
    use pylab to make a histograms of parameter number=param
    if gauss==True, also plot a gaussian with the data's mean and stdev
         (only use the 
    """
    
    try:
        s = MCMC.samples
    except AttributeError:
        s = MCMC

    s1 = s.transpose()[param]
        
    hist = pylab.hist(s1,  bins=nbins, orientation=orientation)

    if gauss:
        mean = s1.mean()
        stdv = stddev(s1)
        print 'mean = %f +- %f' % (mean, stdv)
        smin, smax = min(hist[1]), max(hist[1])
        norm = s1.size * (smax-smin)/nbins  # sum(hist[0])
        ss = numpy.arange(smin, smax, (smax-smin)/5/nbins)
        pylab.plot(ss, norm*GaussPDF(ss, mean, stdv), 'k', hold=True)
        
    return hist

def scatter2d(MCMC, params):
    """
    make a 2d scatterplot of params[0], params[1] from the MCMC.

    TODO: overlay gaussian contours
          color code
          
    """
    try:
        s = MCMC.samples
    except AttributeError:
        s = MCMC

    s1 = s.transpose()[params[0]]
    s2 = s.transpose()[params[1]]

    return pylab.plot(s1,s2, '.')

def histgrid(MCMC, params=None, nbins=30, labels=None):
    """
    make a 2d grid of histograms and scatterplots of the selected parameters
    """

    try:
        s = MCMC.samples
    except AttributeError:
        s = MCMC

    if params is None:
        npar = s.shape[1]
        params=xrange(npar)
    else:
        npar=len(params)

    if labels is None:
        try:
            labels= MCMC.texNames
        except AttributeError:
            pass
        

    nrow = npar
    ncol = npar

    for ipar1, par1 in enumerate(islice(params, 1, npar)):        ### rows 
        for ipar2, par2 in enumerate(islice(params, ipar1+1)):          ### columns
                    
            ax=pylab.subplot(nrow, ncol, npar*ipar1+ipar2+1)
            splot = scatter2d(MCMC, (par2, par1))
                    ### par2, par1 since x-axis along columns
            ax.set_xticklabels([])
            if ipar2 != 0:
                ax.set_yticklabels([])
            elif labels is not None:
                ax.set_ylabel(labels[par1])
                
    for ipar, par in enumerate(params):
        ax=pylab.subplot(nrow, ncol, npar*(npar-1)+ipar+1)
        hist(MCMC, par, nbins=nbins)
        if ipar != 0:
            ax.set_yticklabels([])
        if labels is not None:
            ax.set_xlabel(labels[par])


def convertSampleFile(filename):
    fil = open(filename, 'r')

    lnLike = []; samples = []
    for lin in fil:
        cols = lin.split()
        nsamp = int(cols[0])
        npar = len(cols) - 2
        lnLike.extend([float(cols[1])]*nsamp)
        samples.extend([float(c) for c in cols[2:]]*nsamp)

    fil.close()
    samples = numpy.array(samples)
    nrow = len(samples)/npar
    samples.shape=(nrow, npar)

    return numpy.array(lnLike), samples


def stddev(arr):
    return arr.std()*numpy.sqrt((arr.size-1)/arr.size)


def main(filename, burn = None):
    lnLike, samples = convertSampleFile(filename)

    labels = ['EE1', 'EE2', 'EE3', 'BB1', 'BB2', 'BB3']

    print 'lnLike: ', lnLike.shape, lnLike.dtype
    print 'samples: ', samples.shape, samples.dtype
    histgrid(samples[burn:], labels = labels)
               

if __name__=="__main__":
    main(sys.argv[1])
