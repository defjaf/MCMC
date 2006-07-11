""" deal with 2d beam data (from MAXIPOL) """

from __future__ import division
import sys
import math
import operator
from copy import copy
import os.path
import pickle
import pylab
import gzip
from itertools import izip, repeat

#from numarray import (array, float64, zeros, ones, int32, log, where, ravel, exp,
#                      arange, asarray, sqrt)
#from numarray import concatenate as cat
#import numarray
#import numarray.ma as ma
#from numarray.random_array import uniform

import numpy
from numpy import  (array, float64, zeros, ones, int32, log, where, exp,
                    arange, asarray, sqrt, minimum, maximum, logical_and)
from numpy import concatenate as cat

from numpy.random import uniform
import numpy.core.ma as ma

import MCMC
import Likelihood
from BeamFit import (BeamModel, BeamData, NormalizedBeamModel, 
                     OffsetNormalizedBeamModel, NormalizedBeamLikelihood)
from BeamData import BeamData

import getdist

def readMAXIPOLdataLuis(filename):
    """ read the data and make an instance of class BeamData """

    ia=[]; ja=[]
    i=[]; j=[]; beam=[]; sig=[]; cts=[]
    for line in open(filename, "r"):
        line = line.strip().split()
        i1, j1, b1, s1, c1 = (int(line[0]), int(line[1]), 
           float(line[2]), float(line[3]), int(line[4]))
        ia.append(i1); ja.append(j1)
        if b1 != 0 and s1 !=0:
            ## only keep pixels with data
            i.append(i1); j.append(j1); beam.append(b1)
            sig.append(s1); cts.append(c1)

    beam = asarray(beam, float64)
    sig = asarray(sig, float64)
    ## map i and j (before deletion) onto (-1,1)
    x = array([2*(ii-min(ia))/(max(ia)-min(ia))-1 for ii in i], float64)
    y = array([2*(jj-min(ja))/(max(ja)-min(ja))-1 for jj in j], float64)

    return BeamData(x, y, beam, sig, cts=cts)
        

def readMAXIPOLdataBrad(filename, day=False, sigcut=0.0, ctscut=0, cols=None,
                        nhits=None, neg=False):
    """ 
    read the data and make an instance of class BeamData 
    use Jeff Collins numbers for Day/Night offset if Day==1

    cols (default (2,3)) gives the column number for data and sigma, respectively
    
    - ignore 'counts' column
    """

    if cols is None: cols=(2,3)

    print "Reading data from columns %d-%d" % tuple(cols)
    
    ngood = 0;
    ncut = 0;
    offsets = {'el': 0.295, 'az': 0.05}   ## degrees (brad)
    #offsets = {'el': 0.289, 'az': 0.065}   ## degrees (Jeff)
    az=[]; el=[]; beam=[]; sig=[]; cts=[]
    if filename.endswith('gz'):
        fileobj = gzip.open(filename, "r");
    else:
        fileobj = open(filename, "r");
    for line in fileobj:
        line = line.strip().split()
#        az1, el1, b1, s1, c1 = (
#           float(line[0]), float(line[1]), float(line[2]), float(line[3]), int(line[4]))
        az1, el1, b1, s1, c1= (
           float(line[0]), float(line[1]), float(line[cols[0]]), float(line[cols[1]]),
           long(line[-1]))
        
        if nhits:
            s1 /= sqrt(c1)

        if day:
            az1 += offsets['az']*60.0   ## arcmin
            el1 += offsets['el']*60.0
            
        if s1>sigcut and c1>ctscut:
            ## only keep pixels with good data
            az.append(az1); el.append(el1); beam.append(b1);
            sig.append(s1); cts.append(c1)
            ngood += 1
        else:
            ncut += 1
            

    fileobj.close()
    
    print 'Data read: ncut=%d, ngood=%d' % (ncut, ngood)

    beam = asarray(beam, float64)
    sig = asarray(sig, float64)
    az = asarray(az, float64)
    el = asarray(el, float64)
    cts = asarray(cts, float64)

    if neg is not False and ((neg is None and beam.mean() < 0) or neg):
        print 'negating data'
        beam = -beam

    return BeamData(az, el, beam, sig, cts=cts)


def plot(data):
    """ contour the MAXIPOLBeamData data """
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
        print 'Unknown error: Cannot plot results for this data'

def setup_sampler(dir=None, files=None,  num=None,
                  DayNight=2, LuisBrad=1, useNormalizedBeam=False,
                  cols=(2,3), nhits=None, neg=False, rangeScale=None):
    """
    setup for the sampler for a single beam with data in directory 'dir',
    given by the files in the sequence 'files', or with detector
    'num'.
    
    get the data from column numbers 'cols' from the file
    
    """
    
    numpy.set_printoptions(precision=4, linewidth=150, suppress=True)


#### READ DATA
#    sigcut = 0.2 
#    ctscut = 12  
    sigcut = 0.02  ## minimum sigma for data -- eliminates weird low-error points
    ctscut = 4     ## only keep pixels with this many hits
                   ## nb these also apply to TOI data in which case
                   ## cts is the number of hits in the 'parent' pixel
    
    if dir is None:
        homedir = os.path.expandvars('${HOME}/home')
        if not os.path.exists(homedir):
            homedir = os.path.expandvars('${HOME}')
    
        dir='/'.join((homedir, "cmb/maxipol/beams/"))
    
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

    if LuisBrad == 0:
        data = [ readMAXIPOLdataLuis(dir+fil) for fil in files ]
        xyrange = (-1,1)
    elif LuisBrad == 1:
        if num is None and files is not None: 
            day = [False]*len(files)
            print 'Setting day=', day
                
            data = [
                readMAXIPOLdataBrad(dir+fil, d, sigcut=sigcut, ctscut=ctscut, cols=cols, 
                                    nhits=nhits, neg=neg)
                for fil, d in zip(files, day) ]
        xyrange = data[0]   ### set from the zeroth dataset


#### set likelihood, model

    if useNormalizedBeam:
        #mod = NormalizedBeamModel.NormalizedBeamModel
        #like = NormalizedBeamLikelihood.NormalizedBeamLikelihood(data=data, model=mod)
        #npar = 6

        mod = OffsetNormalizedBeamModel.OffsetNormalizedBeamModel
        like = NormalizedBeamLikelihood.NormalizedBeamLikelihood(data=data, model=mod)
        npar = 9
    else:    
        mod = BeamModel.GaussianBeamModel2D
        like = Likelihood.Likelihood(data=data, model=mod)
        npar = 5

    mod.setxyRange(xyrange, scale=rangeScale)    ## class variables: sets the prior for all instances
    mod.sigMax = 8.0           # 15.0
    mod.sigMin = 3.0

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


def sampleall(nruns=2, nMC=(3000, 100000), useNormalizedBeam=True, irun=0,
              noCorrelations=True, fac=None, doBlock=True):
    """
    run the sampler nruns times for the detectors with both Day and Night data
    """
    
    plotOne = False   ### can't get this to work yet!
    if plotOne:
        pylab.axes()
        pylab.xlim(-100,100)
        pylab.ylim(-100,100)
 
    #dets = [13, 14, 15, 23, 24, 25]   ## LUIS

    if irun==0:
        dets = [12, 13, 14, 15]  # brad, day
        nrow=2; ncol=2
      #  dets = [13, 14] # brad, both
      #  nrow=1; ncol=2
        DayNight=0
    elif irun==1:
        dets = [13, 14, 23, 24, 25, 33, 34, 35, 43, 44, 45 ] # brad, night
        nrow=3; ncol=4
      #  dets = [13, 14] # brad, both
      #  nrow=1; ncol=2
        DayNight = 1
    elif irun==2:    
        dets = [13, 14] # brad, both
        nrow=1; ncol=2
        DayNight=2
    
    reslist = []
    nfig=2
    ntotrun = nruns*nfig
    for run in range(nruns):
        res={}
        for ib, det in enumerate(dets):
            print 'Detector: %d' % det
            fig=pylab.figure(irun*ntotrun+nfig*run)
            if not plotOne:
                ax=fig.add_subplot(nrow, ncol, ib+1)
                ax.cla()

            like, prop_sigmas, start_params = setup_sampler(num=det, DayNight=DayNight,
                                                            useNormalizedBeam=useNormalizedBeam)
            
            res[det] = sample1beam(like, nMC=nMC,  fac=fac,
                                   prop_sigmas=prop_sigmas, start_params=start_params,
                                   noCorrelations=noCorrelations,
                                   doBlock=doBlock)
            if plotOne: 
                pylab.xlim(-100,100)
                pylab.ylim(-100,100)

            fig=pylab.figure(irun*ntotrun+nfig*run+1)
            ax=fig.add_subplot(nrow, ncol, ib+1)
            samples = cat([ s.samples for s in res[det][0] ])
            #samples.transpose()
            ## nb. with numpy, a.transpose doesn't change the array, just gives a new view.
            for var in samples.transpose(): ax.plot(var)
        reslist.append(res)

    return reslist

    
def testTOI(nMC=(3000, 100000), useNormalizedBeam=True,
            noCorrelations=True, fac=None, doBlock=True, cols=None, dets=None,
            mapOnly=False, nhits=None, neg=None, rangeScale=None, 
            closeFigs=False, figName=None, startCols=None, startParams=None):
    """
    run the sampler  for the detectors with TOI data
    """
    
    reslist = []
    nfig=2
    ntotrun = nfig


 #### AHJ: polz only
 ####if dets is None: dets = [13, 14, 15, 23, 24, 25, 33, 34, 35, 43, 44, 45]
 #####   nrow=4; ncol=3
    
#### needed multipliers for data: -, + are strong requirements, 0 is don't care
#                [33, 34, 43, 45]  [13, 14, 15, 23, 24, 25, 35, 44]
#columns 4-5      -    -   -   -     0  0   0+   0   0  0+   0   0-
#columns 6-7     0-   0-   +   -     -  0    0   +   -   -   0   0
    
    alldets = [33, 34, 43, 45, 13, 14, 15, 23, 24, 25, 35, 44]
    neg45   = [-1, -1, -1, -1, +1, +1, +1, +1, +1, +1, +1, -1]
    neg67   = [-1, -1, +1, -1, -1, +1, +1, +1, -1, -1, +1, +1]

    dfac = {}
    dfac[(2,3)] = dict(izip(alldets, repeat(1)))
    dfac[(4,5)] = dict(izip(alldets, neg45))
    dfac[(6,7)] = dict(izip(alldets, neg67))

    neg1 = neg

    if dets is None: dets = [33, 34, 43, 45]

    nrow = ncol = int(math.sqrt(len(dets)))
    if nrow*ncol < len(dets): ncol += 1
    if nrow*ncol < len(dets): nrow += 1

    pref = 'TOI_polarized/b'
    ident = ['_map','_TOI']
    suff = '.txt.gz'
    if mapOnly: 
        ident = [ident[0]] 

    res={}
    for ib, det in enumerate(dets):

        
        res[det] = []
        startres = []
        filebase = [str(det).strip()+id1 for id1 in ident]

        for fb in filebase:
            fil = pref+fb+suff
            figf = fb
            if isinstance(figName, (str, unicode)):
                figf = figName+figf

            print 'Running: ', fil
            fig=pylab.figure(0)
            ax=fig.add_subplot(nrow, ncol, ib+1)
            ax.cla()

            try:
                if startCols is not None:
                    if neg is None: neg1 = (dfac[tuple(startCols)][det]<0)

                    like, prop_sigmas, start_params = setup_sampler(
                        files=[fil],
                        useNormalizedBeam=useNormalizedBeam,
                        cols=startCols, nhits=nhits, neg=neg1,
                        rangeScale=rangeScale)

                    startres, ana = sample1beam(like, nMC=nMC, fac=fac, 
                                           prop_sigmas=prop_sigmas,
                                           start_params=start_params,
                                           noCorrelations=noCorrelations,
                                           doBlock=doBlock)
                                           
                    mod = like.model
                    orig_params = mod.package(ana[0])     ## startres[-1].mean()
                    orig_sigmas = mod.package(3*ana[1])   ## startres[-1].stdev()

                    print ("Start: " + mod.fmtstring) % tuple(mod.unpackage(orig_params))
                    print ("Sigma: " + mod.fmtstring) % tuple(mod.unpackage(orig_sigmas))
                    
                ## need to run this to get the correct likelihood.
                ##   therefore may need to adjust prior ranges

                if neg is None: neg1 = (dfac[tuple(cols)][det]<0)
                like, prop_sigmas, start_params = setup_sampler(
                    files=[fil],
                    useNormalizedBeam=useNormalizedBeam,
                    cols=cols, nhits=nhits, neg=neg1,
                    rangeScale=rangeScale)
                
                if startCols is not None:
                    prop_sigmas = orig_sigmas
                    start_params = orig_params
                    
                if startParams is not None:
                    start_params=startParams
                
                #prev_fig = pylab.gcf()    
                #pylab.figure(4)
                #ll, xx, yy = get_likelihood_grid(like, start_params) 
                ##return ll,xx,yy
                #pylab.figure(prev_fig)

                res[det].append(sample1beam(like, nMC=nMC, prop_sigmas=prop_sigmas,
                                       start_params=start_params, fac=fac, 
                                       noCorrelations=noCorrelations, doBlock=doBlock))
                
                if figName:
                    fig.savefig(figf+str(fig.number).strip()+'.png')

                sys.stdout.flush()
                fig=pylab.figure(1)
                ax=fig.add_subplot(nrow, ncol, ib+1)
                samples = cat([ s.samples for s in res[det][-1][0] ])

                for var in samples.transpose(): ax.plot(var)
                if figName:
                    fig.savefig(figf+str(fig.number).strip()+'.png')

                fig=pylab.figure(2)
                getdist.histgrid(res[det][-1][0][-1])

                if figName:
                    fig.savefig(figf+str(fig.number).strip()+'.png')

            except None:
                print "Unexpected error:", sys.exc_info()[0]
                print "... when running ", fil, fb

            if closeFigs: pylab.close('all')

    pylab.show()
    return res


def saveres(reslist, file=None):
    """
    pickle the results; can't save the actual sampler class since
    you can't pickle all its members
    
    NO LONGER NEEDED; can just pickle the MCMC instances themselves
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

def makereport(reslist, file=sys.stdout, hasTOI=False):
    """
    print out a 'report' of the MCMC run
    """
    for irun, resrun in enumerate(reslist):
        file.write("Run: %d\n" % irun)
        for det, res in resrun.iteritems():
            file.write("%d" % det)

            try:
                if not hasTOI: 
                    val = res[1][0]
                    sig = res[1][1]
                else:
                    val = res[-1][1][0]
                    sig = res[-1][1][1]
                
                for (v, s) in zip(val, sig):
                    file.write("   %f +- %f" % (v, s))
            except:
                print "\n... when running ", irun, det,
                print "Unexpected error:", sys.exc_info()[0]
                    
            file.write("\n")

def plotMod(data, params=None, model=None, hold=False, centeronly=False):
    """ plot the data in MAXIPOLBeamdData.data with params 
    actually, doesn't use the model parameter"""

    x, y, d = regrid(data.x, data.y, data.d)
    #x, y, d = Numeric.asarray(x), Numeric.asarray(y), Numeric.asarray(d)
    if centeronly:   ## doesn't work yet!
        #m2 = where(d < data.sig)   # 1 sig error
        m3 = where(d.filled(0) < max(d.filled(0))*math.exp(-8.0))  ## 4 sig in beam
        newmask = ma.mask_or(d.mask(), m3)
        d = ma.masked_array(d, mask=newmask)
        
    ax = pylab.gca()
    xx, yy = pylab.meshgrid(x,y)
    #ax.pcolor(xx, yy, d, shading='flat', hold='true')
    ax.pcolor(xx, yy, d, shading='flat', hold=hold)
    ## aspect = 'preserve' for square pixels; can't do that with contour however
    if params is not None:
        vals = model(*params).atxy(xx, yy)
        vals.shape = d.shape
        arng=cat(([0.1], arange(5)))
        levs = max(numpy.ravel(vals))*exp(-0.5*arng**2)
        ax.contour(x, y, vals, levs, hold='true', colors='black')
        

class DataSizeError(Exception):
    def __init__(self, *value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        
    
