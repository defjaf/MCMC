from __future__ import division

""" 
    read and process 2d beam data from Planck (based on original from MAXIPOL) 
    
    File info from A Chamballu
    
    White noise levels
    0.000492295 for 143-5 
    0.000404022 for 217-1 
    0.000303889 for 353-1 
    0.000184390 for 545-1

    Here are the names and locations of the maps for the 217-1 bolometer (the more pathologic one): 
     - for the HII region (worst results):
        /space/dmc/m2db/DATA/FPRec_IMG2D_simu4/map_data_217-1_X_Y_HII 
     - for the radio source: 
        /space/dmc/m2db/DATA:FPRec_IMG2D_simu4:map_data_217-1_X_Y_radio 
    with X going from 0 to 19 and Y ( = 1, 2) corresponding to the iteration in the pipeline; 
"""

whitenoise = {143:0.000492295, 217:0.000404022, 353:0.000303889, 545:0.000184390 }


import sys
import math
import os.path
import pickle
import pylab
import gzip
from itertools import izip, repeat

import numpy
from numpy import (array, float64, zeros, ones, int32, log, where, exp,
                   arange, asarray, sqrt, minimum, maximum, logical_and, empty)
from numpy import concatenate as cat

from BeamData import BeamData

from ProcessBeamData import setup_sampler, sample1beam, plotter

import getdist

try:
    import piolib as pio
except:
    print "Will probably fail, piolib not available"

#    sigcut = 0.2 
#    ctscut = 12  
sigcut = 0.02  ## minimum sigma for data -- eliminates weird low-error points
ctscut = 4     ## only keep pixels with this many hits
               ## nb these also apply to TOI data in which case
               ## cts is the number of hits in the 'parent' pixel

def getobjname(db=None, grp=None, det=None, MC=None, iter=None, obj=None):
    if db is None: 
        db = "/space/dmc/m2db/DATA"
    if grp is None: 
        grp = "FPRec_IMG2D_simu4"
        
    base = db + '/' + grp + '/' + 'map_' 
    suff = '_' + str(MC).strip() + '_' + str(iter).strip()+'_' 
    return [base + kind + det + suff + obj for kind in ['data_', 'hit_']]


##MAXI   split from original setup_sampler
def read_data_Planck(files=None, sigcut=0.0, ctscut=0, nhits=None, neg=False, noise=None, **kwargs):
    """
    
    read piolib img2d data. For each detector, read a pair (data, nhit).
    
    if files is present, it is a (2,N) sequence where each entry is (data, nhit)
    
    otherwise, call getobjname(**kwargs) to get the files -- do this, or just explicitly call???
    """
    
    
    if files is None:
        files = getobjname(**kwargs)

    if noise is not None:
        sigma_white = noise
    else:
        try:
            det = kwargs['det']
            sigma_white = sqrt(whitenoise[int(det[:3])])  ## this is a cheat to convert the det to an integer
        except KeyError:   ## should catch appropriate exception
            print "Can't get white noise; using 1"
            sigma_white = 1.0

    data = []
    for iset in files:
        img = pio.ReadIMG2DObject(files[0], "PIODOUBLE", "")
        hit = pio.ReadIMG2DObject(files[1], "PIODOUBLE", "") ### double or int???
        sig = sigma_white/sqrt(hit)
        ## assume square, although really should read keywords for size
        npix = sqrt(img.size)
        if npix*npix != img.size:
            print "Image size problem: size=%d, npix=%f" % (img.size, npix)
        
        ### need to generate (x,y) position arrays. Just use integers?
        
        x = empty((npix, npix),dtype=float64)
        y = empty((npix, npix),dtype=float64)
         
        x[:] = arange(npix)
        y[:] = arange(npix)
        
        y = y.T
        
        x = asarray(x, float64).flatten()
        y = asarray(y, float64).flatten()
        img = asarray(img, float64).flatten()
        sig = asarray(sig, float64).flatten()
        hit = asarray(hit, float64).flatten()
        
        data.append(BeamData(x, y, img, sig, cts = hit))
        
    xyrange = data[0]   ### set from the zeroth dataset

    return data, xyrange


## rewritten for Planck -- just the very inner loop is generic
def sampleall(nruns=2, nMC=(3000, 100000), useNormalizedBeam=True, irun=0,
              noCorrelations=True, fac=None, doBlock=True, nhits=None):
    """
    run the sampler nruns times 
    """
    
    dets = ["217-1"]
    objs = ["HII", "Radio"]
    iters = range(1,3)
    MCs = range(20)
    
    
    
    reslist = []
    nfig=2
    ntotrun = nruns*nfig
    for run in range(nruns):
        res={}
        ib = 0
        for det in dets:
            for obj in objs:
                for it in iters:
                    for MC in MCs:
                        ib += 1
                        print 'Detector: %d, obj: %s, iter: %d, MC: %d' % (det, obj, it, MC)
                        fig=pylab.figure(irun*ntotrun+nfig*run)
                        ax=fig.add_subplot(nrow, ncol, ib+1)
                        ax.cla()

                        ## need to explicitly read the data here, now -- how to make generic?
                        data, xyrange = read_data_Planck(det=det, MC=MC, iter=it, obj=obj, 
                                                         sigcut=sigcut, ctscut=ctscut)
                        like, prop_sigmas, start_params = setup_sampler(data, xyrange,
                                                                        useNormalizedBeam=useNormalizedBeam,
                                                                        sigminmax=(2,20))
            
                        res[det] = sample1beam(like, nMC=nMC,  fac=fac,
                                               prop_sigmas=prop_sigmas, start_params=start_params,
                                               noCorrelations=noCorrelations,
                                               doBlock=doBlock)

                        fig=pylab.figure(irun*ntotrun+nfig*run+1)
                        ax=fig.add_subplot(nrow, ncol, ib+1)
                        samples = cat([ s.samples for s in res[det][0] ])
                        for var in samples.transpose(): ax.plot(var)
        reslist.append(res)

    return reslist


#  corrected for Planck ##MAXI
def testPlanck(nMC=(3000, 100000), useNormalizedBeam=True,
            noCorrelations=True, fac=None, doBlock=True, dets=None,
            nhits=None, rangeScale=None, 
            closeFigs=False, figName=None, startParams=None):
    """
    run the sampler  for the detectors with TOI data
    """
    
    reslist = []
    nfig=2
    ntotrun = nfig

    
    if dets is None: dets = ["217-1"]
    objs = ["HII", "Radio"]
    iters = range(1,3)
    MCs = range(20)


    ntot = len(dets) * len(objs) * len(iters) * len(MCs)
    nrow = ncol = int(math.sqrt(ntot))
    if nrow*ncol < ntot: ncol += 1
    if nrow*ncol < ntot: nrow += 1


    res={}
    ib = 0
    for det in dets:
        for obj in objs:
            for it in iters:
                for MC in MCs:
                    print 'Detector: %s, obj: %s, iter: %d, MC: %d' % (det, obj, it, MC)
                    
                    res[ib] = []
                    startres = []

                    fig=pylab.figure(0)
                    ax=fig.add_subplot(nrow, ncol, ib+1)
                    ax.cla()
        
                    try:                
                        ## need to run this to get the correct likelihood.
                        ##   therefore may need to adjust prior ranges

                        data, xyrange = read_data_Planck(det=det, MC=MC, iter=it, obj=obj, 
                                                         sigcut=sigcut, ctscut=ctscut)
                        like, prop_sigmas, start_params = setup_sampler(
                            data, xyrange,
                            useNormalizedBeam=useNormalizedBeam,sigminmax=(2,20),
                            rangeScale=rangeScale)
                            
                        if startParams is not None:
                            start_params=startParams
            

                        res[ib].append(sample1beam(like, nMC=nMC, prop_sigmas=prop_sigmas,
                                               start_params=start_params, fac=fac, 
                                               noCorrelations=noCorrelations, doBlock=doBlock))
            
                        if figName:
                            fig.savefig(figf+str(fig.number).strip()+'.png')

                        sys.stdout.flush()
                        fig=pylab.figure(1)
                        ax=fig.add_subplot(nrow, ncol, ib+1)
                        samples = cat([ s.samples for s in res[ib][-1][0] ])

                        for var in samples.transpose(): ax.plot(var)
                        if figName:
                            fig.savefig(figf+str(fig.number).strip()+'.png')

                        fig=pylab.figure(2)
                        getdist.histgrid(res[ib][-1][0][-1])

                        if figName:
                            fig.savefig(figf+str(fig.number).strip()+'.png')

                    except None:
                        print "Unexpected error:", sys.exc_info()[0]
                        print "... when running Detector: %d, obj: %s, iter: %d, MC: %d" % (det, obj, it, MC)
                        
                    ib += 1
                    

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


def makereport(reslist, file=sys.stdout, hasTOI=False, hasRuns=False):
    """
    print out a 'report' of the MCMC run
    """

    if not hasRuns: reslist=[reslist]
    
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

    
