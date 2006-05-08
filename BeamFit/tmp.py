""" deal with 2d beam data (from MAXIPOL) """

from __future__ import division
import sys
import math
import operator
from copy import copy
import os.path
import pickle
import pylab

from numarray import array, Float64, zeros, ones, Int, log
from numarray import concatenate as cat
import numarray
import numarray.ma as ma
from numarray.random_array import uniform

#from scipy.base import array, Float64, zeros, ones, Int, log, empty
#from scipy.base import concatenate as cat
#import scipy.base
#import scipy.base.ma as ma
#from scipy.random_array import uniform

import MCMC
import Likelihood
from BeamFit import BeamModel, BeamData, NormalizedBeamModel, NormalizedBeamLikelihood
from BeamData import BeamData


def readMAXIPOLdataLuis(filename):
    """ read the data and make an instance of class BeamData """

    ia=[]; ja=[]
    i=[]; j=[]; beam=[]; sig=[]; cts=[]
    for line in open(filename, "r"):
        line = line.strip().split()
        i1, j1, b1, s1, c1 = \
           int(line[0]), int(line[1]), float(line[2]), float(line[3]), int(line[4])
        ia.append(i1); ja.append(j1)
        if b1 != 0 and s1 !=0:
            ## only keep pixels with data
            i.append(i1); j.append(j1); beam.append(b1); sig.append(s1); cts.append(c1)

    beam = array(beam, Float64)
    sig = array(sig, Float64)
    ## map i and j (before deletion) onto (-1,1)
    x = array([2*(ii-min(ia))/(max(ia)-min(ia))-1 for ii in i], Float64)
    y = array([2*(jj-min(ja))/(max(ja)-min(ja))-1 for jj in j], Float64)

    return BeamData(x, y, beam, sig)
        

def readMAXIPOLdataBrad(filename, day=False, sigcut=0.0):
    """ 
    read the data and make an instance of class BeamData 
    use Jeff Collins numbers for Day/Night offset if Day==1
    """

    offsets = {'el': 0.289, 'az': 0.065}   ## degrees
    az=[]; el=[]; beam=[]; sig=[]; cts=[]
    for line in open(filename, "r"):
        line = line.strip().split()
        az1, el1, b1, s1, c1 = \
           float(line[0]), float(line[1]), float(line[2]), float(line[3]), int(line[4])
        if day:
            az1 -= offsets['az']*60.0   ## arcmin
            el1 -= offsets['el']*60.0
        if s1>sigcut:
            ## only keep pixels with data
            az.append(az1); el.append(el1); beam.append(b1); sig.append(s1); cts.append(c1)

    beam = array(beam, Float64)
    sig = array(sig, Float64)
    az = array(az, Float64)
    el = array(el, Float64)

    return BeamData(az, el, beam, sig)


def plot(data):
    """ contour the MAXIPOLBeamData data """
    x, y, d = regrid(data.x, data.y, data.d)
    pylab.imshow(ma.filled(d,0), extent=[min(x), max(x), min(y), max(y)],
                 interpolation='nearest', origin='lower', aspect='free')
    return pylab.contour(x, y, ma.log(d))    


def regrid(x, y, data, fill=0):
    """
    assuming that x, y are a subset of uniformly-gridded points,
    return a fully-populated 2-d array of data and 1-d arrays of x, y

    fill unavailable data with fill

    [can work with a sequence of arrays in data[] ** not yet** ]
    """

    xx = grid1d(x)
    yy = grid1d(y)

    d = fill + zeros(shape=(len(yy), len(xx)), type=Float64)
    mask = ones(shape=(len(yy), len(xx)), type=Int)
    for x1, y1, d1 in zip(x, y, data):
        i = int(round((x1-xx[0])/(xx[1]-xx[0])))
        j = int(round((y1-yy[0])/(yy[1]-yy[0])))
        d[j, i] = d1
        mask[j,i] = 0

        dat=ma.array(d, mask=mask)

        ## rewrite to do a sequence of d[:, i, j]???

    return  xx, yy, dat
    

def grid1d(x, dx=0, nx=0):
    """
    make an evenly-spaced 1d grid of x values out of the subsequence of given x
    """
    
    u = []
    for e in x:
        if e not in u:
            u += [e]

    u.sort()
    nu = len(u)
    minx = min(u)
    maxx = max(u)
    diff = [u[i+1]-u[i] for i in range(nu-1)]
    dx = min(diff)
    n = int(round((max(u)-min(u))/dx)+1)

    return [minx + dx*i for i in range(n)] 

def plotter(sampler):
    """
    make a plot from the results of the sampler
    (could probably make this more generic if plotMod is a member fn of the model)
    for now, only plot the *first* dataset (since for this case that's the 'good' data)
    """
    
    mod = sampler.like.model
    data = sampler.like.data

    ntot = len(sampler.samples)
    stride = max(1,ntot//sampler.naccept)
    ana = MCMC.chain_analyze(sampler.samples[(ntot//5)::stride,:])
    vals = sampler.like.model.package(ana[0])
    sigs = sampler.like.model.package(ana[1])
    print vals
    print sigs
    print ana[2]
        
    pylab.cla()
    #for d in data: plotMod(d, vals, mod)
    plotMod(data[0], vals, mod)

def sample1beam(dir=None, files=None, nMC=(1000,1000), num=None,
                DayNight=2, LuisBrad=1, plotRes=None, useNormalizedBeam=False):
    """
    run the sampler for a single beam with data in directory 'dir',
    given by the files in the sequence 'files', or with detector
    'num'.
    
    do a sequence of MCMC runs with the number of samples given in nMC
    (start each one at the mean of the previous using the previous
    covariance, scaled appropriately, for a gaussian sampling
    distribution)

    if a single run fails (because no new samples are accepted, or the
    variance is zero, continue from there with the same number of
    samples to be added)

    return the full MCMC sample class as well as the summary
    statistics for the last run

    set DayNight=0,1,2 for Day, Night, Both
    """
    
    sigcut = 0.02  ## minimum sigma for data -- eliminates weird low-error points
    
    if dir is None:
        dir=os.path.expandvars("${HOME}/cmb/maxipol/beams/")
    if num is not None:
        if LuisBrad == 0: # Luis names
            filen = "binnedb"+str(num).strip()+".txt"
            if DayNight == 0:
                files=["Day/"+filen]
            elif DayNight == 1:
                files=["Night/"+filen]
            elif DayNight == 2:
                files=["Day/"+filen, "Night/"+filen]
        elif LuisBrad == 1: # Brad names
            filen = "brad/b"+str(num).strip()
            suff = 'time_jupiter_maps.txt'
            dnstring = ["day", "night"]
            ## array day is True for daytime, False for nighttime (see read*Brad above)
            if DayNight == 0 or DayNight == 1:
                files=[filen+dnstring[DayNight]+suff]
                day = [not DayNight]
            elif DayNight == 2:
                files=[filen+dnstring[id]+suff for id in [0,1]]
                day = [True, False]
            
    if files is None:
        files=["Day/binnedb25.txt", "Night/binnedb25.txt"]
        day = [True, False]

    if operator.isNumberType(nMC):
        nMC = (nMC, nMC)
    if LuisBrad == 0:
        data = [ readMAXIPOLdataLuis(dir+fil) for fil in files ]
        xyrange = (-1,1)
    elif LuisBrad == 1:
        data = [ readMAXIPOLdataBrad(dir+fil, d, sigcut=sigcut) for fil, d in zip(files, day) ]
        xyrange = data[0]   ### set from the zeroth dataset


    if useNormalizedBeam:
        mod = NormalizedBeamModel.NormalizedBeamModel
        like = NormalizedBeamLikelihood.NormalizedBeamLikelihood(data=data, model=mod)
        npar = 6
    else:    
        mod = BeamModel.GaussianBeamModel2D
        like = Likelihood.Likelihood(data=data, model=mod)
        npar = 5
        
    mod.setxyRange(xyrange)    ## class variables: sets the prior for all instances
    mod.setsigMax(xyrange)
    

    dx = (mod.centerMin[0],mod.centerMax[0])
    dy = (mod.centerMin[1],mod.centerMax[1])
    delx = dx[1]-dx[0]
    dely = dy[1]-dy[0]
    
    prop_sigmas = ( (delx/3, dely/3 ), (delx/10, dely/10 ), 0.3)

 #   start_params = ( (uniform(*dx), uniform(*dy)), 
 #                    (uniform(0,delx)/5, uniform(0,dely)/5),
 #                    uniform(0,math.pi/2) )
    start_params = ( ((dx[0]+dx[1])/2, (dy[0]+dy[1])/2), 
                     (delx/5, dely/5), 0) 

    if useNormalizedBeam:
        prop_sigmas += max(data[0].d)/4,   ## nb commas for tuples; "+=" appends
        start_params += max(data[0].d)/2,
        
    if plotRes is None:
        return MCMC.sampler(like, nMC, prop_sigmas, start_params, plotter=plotter)
    else:
        ## also make it work on the pickled data
        ana = plotRes[-1][1]
        vals = like.model.package(ana[0])
        sigs = like.model.package(ana[1])
        
        pylab.cla()
        for d in data: plotMod(d, vals, mod)


def sampleall(nruns=2, nMC=(3000, 100000), useNormalizedBeam=False, irun=0):
    """
    run the sampler nruns times for the detectors with both Day and Night data
    """

    #dets = [13, 14, 15, 23, 24, 25]   ## LUIS

    if irun==0:
        dets = [12, 13, 14, 15, 33]  # brad, day
        nrow=2; ncol=3
        DayNight=0
    elif irun==1:
        dets = [13, 14, 15, 23, 24, 25, 34, 35, 43, 44, 45 ] # brad, night
        nrow=3; ncol=3
        DayNight = 1
    elif irun==2:    
        dets = [13, 14, 15] # brad, both
        nrow=2; ncol=2
        DayNight=2
    
    reslist = []
    nfig=2
    ntotrun = nruns*nfig
    for run in range(nruns):
        res={}
        for ib, det in enumerate(dets):
            print 'Detector: %d' % det
            pylab.figure(irun*ntotrun+nfig*run)
            pylab.subplot(nrow, ncol, ib+1)
            pylab.cla()
            res[det] = sample1beam(num=det, nMC=nMC, DayNight=DayNight,
                       useNormalizedBeam=useNormalizedBeam)
            pylab.figure(irun*ntotrun+nfig*run+1)
            pylab.subplot(nrow, ncol, ib+1)
            samples = cat([ s.samples for s in res[det][0] ])
            samples.transpose()
            for var in samples: pylab.plot(var)
        reslist.append(res)

    return reslist

def saveres(reslist, file=None):
    """
    pickle the results; can't save the actual sampler class since
    you can't pickle all its members
    """
    newres=[]
    for resrun in reslist:
        for det, res in resrun.iteritems():
            resrun[det]=(res[0][-1].samples, res[1])  ## save the 'ana' element & samples
            ## should save in a form closer to the original???
        newres.append(resrun)
    if file is not None:
        fp = open(file, 'w')
        pickle.dump(newres, fp)
        fp.close()
    return newres

def makereport(reslist, file=sys.stdout):
    """
    print out a 'report' of the MCMC run
    """
    for irun, resrun in enumerate(reslist):
        file.write("Run: %d\n" % irun)
        for det, res in resrun.iteritems():
            file.write("%d" % det)
            
            val = res[1][0]
            sig = res[1][1]
            for (v, s) in zip(val, sig):
                file.write("   %f +- %f" % (v, s))
            file.write("\n")

def samplef(nruns=2, nMC=None):
    """ 
    run the sampler nruns times for the detectors with an explicit list of files
    >>>>> not used
    """

  #  dets = [13, 14, 15, 23, 24, 25]
    files=[["Day/binnedb14.txt"]]
    dets =[14]
    nrow = 2; ncol = 3
    if nMC is None: nMC = (3000, 100000) ##  (500, 1000)
    reslist = []
    nfig=2
    for run in range(nruns):
        res={}
        for ib, (det, fil) in enumerate(zip(dets, files)):
            pylab.figure(nfig*run)
            pylab.subplot(nrow, ncol, ib+1)
            pylab.cla()
            res[det] = sample1beam(files=fil, nMC=nMC)
            pylab.figure(nfig*run+1)
            pylab.subplot(nrow, ncol, ib+1)
            samples = cat([ s.samples for s in res[det][0] ])
            samples.transpose()
            for var in samples: pylab.plot(var)
        reslist.append(res)

    return reslist

def plotMod(data, params=None, model=None, hold=False):
    """ plot the data in MAXIPOLBeamdData.data with params 
    actually, doesn't use the model parameter"""
    x, y, d = regrid(data.x, data.y, data.d)
    ### make full 2d x, y arrays (there's probably a more clever way to do this!)
    ij = 0
    xx = array(shape=len(x)*len(y), type=Float64)
    yy = array(shape=len(x)*len(y), type=Float64)
    for j in range(len(y)):
        for i in range(len(x)):
            xx[ij] = x[i]
            yy[ij] = y[j]
            ij += 1
    pylab.imshow(ma.filled(d,0), extent=[min(x), max(x), min(y), max(y)],
                     interpolation='nearest', origin='lower', aspect='free', hold=hold)
    ## aspect = 'preserve' for square pixels; can't do that with contour however
    #pylab.contour(x, y, ma.log(d))
    if params is not None:
        vals = model(*params).atxy(xx, yy)
        vals.shape = d.shape
        pylab.contour(x, y, vals)
        
