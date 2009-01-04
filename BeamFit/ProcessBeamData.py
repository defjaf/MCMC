""" deal with 2d beam data (originally from MAXIPOL) """

## refactored from MAXIPOLBeamData???
## see ##MAXI below for candidates for latter

from __future__ import division
import sys
import math
import operator
import pylab

import numpy
from numpy import  (array, float64, zeros, ones, int32, log, where, exp,
                    arange, asarray, sqrt, minimum, maximum, logical_and, empty)
from numpy import concatenate as cat

try:
    import numpy.ma as ma
except ImportError:
    import numpy.core.ma as ma
    
import MCMC
import Likelihood
from BeamFit import (BeamModel, NormalizedBeamModel, 
                     OffsetNormalizedBeamModel, NormalizedBeamLikelihood)

import getdist



def plot(data):
    """ contour the Beam data """
    x, y, d = regrid(data.x, data.y, data.d)
    return pylab.gca().imshow(ma.filled(d,0), extent=[min(x), max(x), min(y), max(y)],
                 interpolation='nearest', origin='lower', aspect='free')
    #return pylab.contour(x, y, ma.log(d))    


def regrid(x, y, data, fill=0):
    """
    assuming that x, y are a subset of uniformly-gridded points,
    return a fully-populated 2-d array of data and 1-d arrays of x, y

    fill unavailable data with fill

    [can work with a sequence of arrays in data[] ** not yet** ]
    
    TODO: more efficient rebinning that doesn't require shape=(lxx,lyy)
    """

    xx = grid1d(x)
    yy = grid1d(y)
    
    lxx = len(xx); lyy = len(yy)

    if 8*lxx*lyy/1024/1024 > 256:
        raise DataSizeError        

    d = fill + zeros(shape=(lyy, lxx), dtype=float64)
    mask = ones(shape=(lyy, lxx), dtype=int32)
    for x1, y1, d1 in zip(x, y, data):
        i = int(round((x1-xx[0])/(xx[1]-xx[0])))
        j = int(round((y1-yy[0])/(yy[1]-yy[0])))
        d[j, i] = d1
        mask[j,i] = 0

        dat=ma.array(d, mask=mask)

        ## rewrite to do a sequence of d[:, i, j]???

    return asarray(xx), asarray(yy), dat
    

def grid1d(x, dx=0, nx=0):
    """
    make an evenly-spaced 1d grid of x values out of the subsequence of given x
    """
            
    u = list(set(x))   ## remove duplicates
    u.sort()
    nu = len(u)
    minx = min(u)
    maxx = max(u)
    diff = [u[i+1]-u[i] for i in range(nu-1)]
    dx = min(diff)
    n = int(round((max(u)-min(u))/dx)+1)

    return [minx + dx*i for i in range(n)] ### use linspace?

def plotter(sampler):
    """
    make a plot from the results of the sampler
    (could probably make this more generic if plotMod is a member fn of the model)
    for now, only plot the *first* dataset (since for this case that's the 'good' data)
    """
    
    mod = sampler.like.model
    data = sampler.like.data

    params = where(sampler.prop.sigmas>0)[0] # (where returns a tuple!)

    ntot = len(sampler.samples)
    stride = 1 #max(1,ntot//sampler.naccept)
    ana = MCMC.chain_analyze(sampler.samples[(ntot//5)::stride,:], params=params)
    vals = sampler.like.model.package(ana[0])
    sigs = sampler.like.model.package(ana[1])
    print sampler.fmtstring % tuple(ana[0])
    print sampler.fmtstring % tuple(ana[1])
    print ana[2]
        
    pylab.cla()
    #for d in data: plotMod(d, vals, mod)
    try:
        plotMod(data[0], vals, mod)
    except DataSizeError:
        print 'Data too big: Cannot plot results for this data'
    except:
        print 'Unknown error: Cannot plot results for this data:'
        print sys.exc_info()[0]



def setup_sampler(data, xyrange, useNormalizedBeam=False ,
                    rangeScale=None, sigminmax=(3,8)):
    """
    setup the sampler using data over the range xyrange
    """

    numpy.set_printoptions(precision=4, linewidth=150, suppress=True)

#### set likelihood, model

    if useNormalizedBeam:
        #mod = NormalizedBeamModel.NormalizedBeamModel
        #like = NormalizedBeamLikelihood.NormalizedBeamLikelihood(data=data, model=mod)
        #npar = 6

        mod = OffsetNormalizedBeamModel.OffsetNormalizedBeamModel
        like = NormalizedBeamLikelihood.NormalizedBeamLikelihood(data=data, model=mod)
        npar = 9
        use_xy = OffsetNormalizedBeamModel.use_xy
    else:    
        mod = BeamModel.GaussianBeamModel2D
        like = Likelihood.Likelihood(data=data, model=mod)
        npar = 5
        use_xy = False

    mod.setxyRange(xyrange, scale=rangeScale)    ## class variables: sets the prior for all instances
    mod.sigMin, mod.sigMax=sigminmax
    print 'setting sigMin, sigMax=', mod.sigMin, mod.sigMax

    dx = (mod.centerMin[0],mod.centerMax[0])
    dy = (mod.centerMin[1],mod.centerMax[1])
    delx = dx[1]-dx[0]
    dely = dy[1]-dy[0]

    print 'center min, max=', mod.centerMin, mod.centerMax


### set initial width and location of proposal

    need_prop_sig = need_start_params = True
#    need_prop_sig = need_start_params = False
#    if prop_sigmas is None:
#        need_prop_sig = True
#    if start_params is None:
#        need_prop_sig = True

    if need_prop_sig:
        # prop_sigmas = ( (delx/3, dely/3 ), (delx/5, dely/5 ), 0.6)
        prop_sigmas = ( (delx/10, dely/10), (delx/10, dely/10), 0.6)
        if use_xy:
            prop_sigmas[-1] = 0.2


 #   start_params = ( (uniform(*dx), uniform(*dy)), 
 #                    (uniform(0,delx)/5, uniform(0,dely)/5),
 #                    uniform(0,math.pi/2) )

    print "averaging over all data for starting point"
    #stats = data[0].stats(sigcut=None, do_abs=False)
    stats = data[0].stats(sigcut=0.0)

    #start_params = ( ((dx[0]+dx[1])/2, (dy[0]+dy[1])/2), 
    #                 (delx/5, dely/5), 0) 

    if need_start_params:
        startsigs = minimum([sqrt(stats[2]), sqrt(stats[3])], mod.sigMax)
        startsigs = maximum(startsigs, mod.sigMin)
        start_params = ((stats[0], stats[1]), startsigs, 0)

    if useNormalizedBeam:
        ## amplitude
        if need_prop_sig:           ## nb commas for tuples; "+=" appends

            # no need to be this fancy, really, just make sure that sigma is large 
            #  enough to "jump the barrier"!
            x = data[0].x-stats[0]; y=data[0].y-stats[1]; d=data[0].d
            dctr = d[logical_and(abs(x)<startsigs[0], abs(y)<startsigs[1])]

            prop_sigmas += 2*(max(dctr)-min(dctr)),   # used to be /10 or something -- bad!
            prop_sigmas += (0.3, 0, 0),   # offset, gradient [cos(theta), phi]
            #   prop_sigmas += (0, 0, 0),

        if need_start_params:
            start_params += (max(data[0].d)-min(data[0].d))/2,
            start_params += (0, 1, 0),  # offset, gradient [cos(theta), phi]


        ### remove the param blocks corresponding to the sig=0 parameters
        ### can't do this with -= etc
        ###            since those chase 'singleton' class variable each time!

        mod.paramBlocks = [0,1,2,3,4,5,6]
        mod.nBlock = 7
#        mod.paramBlocks = [0,1,2,3,4,5,6,7,8]
#        mod.nBlock = 9
   #     mod.paramBlocks = [0,1,2,3,4,5]
   #     mod.nBlock = 6

    print ("Starting point:  " + mod.fmtstring) % tuple(mod.unpackage(start_params))
    print ("Starting sigmas: " + mod.fmtstring) % tuple(mod.unpackage(prop_sigmas))

    ### needs to
    return like, prop_sigmas, start_params


def get_likelihood_grid(like, params):
    """calculate the likelihood in an x,y grid with other params set as in params
    """
    xx = numpy.linspace(like.data[0].x.min(), like.data[0].x.max())
    yy = numpy.linspace(like.data[0].y.min(), like.data[0].y.max())
    xg,yg = pylab.meshgrid(xx,yy)
    
    def mylike(x1,y1):
        myparams = [(x1,y1)]
        myparams.extend(params[1:])
        return like(tuple(myparams))
    
    ll = numpy.array([mylike(x,y) for x in xx for y in yy])
    ll.shape=(len(xx),len(yy))
    ax = pylab.gca()
    #ax.pcolor(xx, yy, d, shading='flat', hold='true')
    ax.pcolor(xg, yg, numpy.transpose(ll), shading='flat')
    pylab.show()
    return xg, yg, ll
                

def sample1beam(like,  prop_sigmas, start_params, nMC=(1000,1000), 
                fac=None, plotRes=None, noCorrelations=False, doBlock=False):
    """
    return the full MCMC sample class as well as the summary
    statistics for the last run
    
    do a sequence of MCMC runs with the number of samples given in nMC
    (start each one at the mean of the previous using the previous
    covariance, scaled appropriately, for a gaussian sampling
    distribution)

    if a single run fails (because no new samples are accepted, or the
    variance is zero, continue from there with the same number of
    samples to be added)
    
    return the likelihood object, the starting parameters, 
    and the starting variances for the proposal

    set DayNight=0,1,2 for Day, Night, Both
    """
    
    #plotter = None
    #print "not using plotter for now"
    if operator.isNumberType(nMC):
         nMC = (nMC, nMC)
             
    if plotRes is None:
        return MCMC.sampler(like, nMC, prop_sigmas, start_params, 
                            plotter=plotter, fac=fac, 
                            noCorrelations=noCorrelations, doBlock=doBlock)
    else:
        ## also make it work on the pickled data
        ana = plotRes[-1][1]
        vals = like.model.package(ana[0])
        sigs = like.model.package(ana[1])
        
        pylab.gca().cla()
        for d in data: plotMod(d, vals, mod)


def plotMod(data, params=None, model=None, hold=False):
    """ plot the data with params 
    actually, doesn't use the model parameter"""
    x, y, d = regrid(data.x, data.y, data.d)
    ### make full 2d x, y arrays (there's probably a more clever way to do this!)
    ij = 0
    xx = empty(len(x)*len(y), dtype=float64)
    yy = empty(len(x)*len(y), dtype=float64)
    for j in range(len(y)):
        for i in range(len(x)):
            xx[ij] = x[i]
            yy[ij] = y[j]
            ij += 1
    pylab.imshow(ma.filled(d,0), extent=[min(x), max(x), min(y), max(y)],
                 interpolation='nearest', origin='lower', aspect='auto', hold=hold)
    ## aspect = 'preserve' for square pixels; can't do that with contour however
    #pylab.contour(x, y, ma.log(d))
    if params is not None:
        vals = model(*params).atxy(xx, yy)
        vals.shape = d.shape
        pylab.contour(x, y, vals)




class DataSizeError(Exception):
    def __init__(self, *value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        
    
