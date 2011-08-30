from __future__ import division

""" 
    read and process TOI data from Planck (based on original from MAXIPOL
       modified for Planck IMG2D data). This version uses flat files from A Chamballu 
    
    need white noise level from data RMS. (Not ideal)
    
    mostly works. todo:
    
       check paramblock functionality. better to do all in one block after 1st run?
       better proposal width?
    

"""

import sys
import math
import os.path
import matplotlib.pyplot as plt
import itertools as itt

import numpy as np
from numpy import (array, float64, zeros, ones, int32, log, where, exp, linspace,
                   arange, asarray, sqrt, minimum, maximum, logical_and, empty)
from numpy import concatenate as cat

from BeamData import BeamData
from ProcessBeamData import setup_sampler, sample1beam, plotter
import getdist

sigminmax=[0,10]

               
def read_data_Planck_TOI(files=None, sigma=None):
    """
    read TOI data. 
    """
    
    data = []
    for fset in files:
        x, y, img, model = np.loadtxt(fset, unpack=True)
    
        if sigma is None:
            sigma_white = np.std(img)
            print "%s: setting sigma=%f" % (fset,sigma_white)
        else:
            sigma_white = sigma

        data.append(BeamData(x, y, img, sigma_white))
        
    xyrange = data[0]   ### set from the zeroth dataset

    return data, xyrange


### not sure what the difference between "sampleall" and "testPlanck". 
###   the latter explicitly makes a plot...

## rewritten for Planck -- just the very inner loop is generic
def sampleall(nruns=2, nMC=(3000, 100000), useNormalizedBeam=True, irun=0,
              noCorrelations=True, fac=None, doBlock=True, nhits=None, rotateParams=False):
    """
    run the sampler nruns times 
    """
    
    iters = range(1,3)
    MCs = range(20)
    
    fdir = os.path.expanduser("~/FPtesting/Beams/TOIs/")
    files = ["model_10_143_5.dat", "model_84_217_4.dat"]
    
    
    reslist = []
    nfig=2
    ntotrun = nruns*nfig
    for run in range(nruns):
        res=[]
        ib = 0
        for f in files:
            for it in iters:
                for MC in MCs:
                    ib += 1
                    print 'file: %s, iter: %d, MC: %d' % (f, it, MC)
                    fig=plt.figure(irun*ntotrun+nfig*run)
                    ax=fig.add_subplot(nrow, ncol, ib+1, adjustable='box',aspect='equal')
                    ax.cla()

                    data, xyrange = read_data_Planck_TOI(f[dir+f])
                    
                    like, prop_sigmas, start_params = setup_sampler(
                        data, xyrange, useNormalizedBeam=useNormalizedBeam,
                                      sigminmax=sigminmax)
        
                    res.append(sample1beam(like, nMC=nMC,  fac=fac,
                                           prop_sigmas=prop_sigmas, start_params=start_params,
                                           noCorrelations=noCorrelations, rotateParams=rotateParams,
                                           doBlock=doBlock))

                    fig=plt.figure(irun*ntotrun+nfig*run+1)
                    ax=fig.add_subplot(nrow, ncol, ib+1)
                    samples = cat([ s.samples for s in res[-1][0] ])
                    for var in samples.transpose(): ax.plot(var)
        reslist.append(res)

    return reslist


#  corrected for Planck ##MAXI
def testPlanck(nMC=(3000, 10000), useNormalizedBeam=True,
            noCorrelations=True, fac=None, doBlock=True, 
            MCs=None, rotateParams=False,
            nhits=None, rangeScale=None, 
            closeFigs=False, figName=None, startParams=None):
    """
    run the sampler  for the TOI data
    """
    
    reslist = []
    nfig=2

    if MCs is None: MCs = [1]

    fdir = os.path.expanduser("~/FPtesting/Beams/TOIs/")
    files = ["model_10_143_5.dat", "model_84_217_4.dat"]

    ntot = len(files) * len(MCs)
    nrow = ncol = int(math.sqrt(ntot))
    if nrow*ncol < ntot: ncol += 1
    if nrow*ncol < ntot: nrow += 1

    res={}
    for ib, (f, MC) in enumerate(itt.product(files, MCs)):
        print 'File: %s, MC: %d' % (f, MC)
        figf = '_'.join([figName, f, str(MC).strip()])
        
        res[ib] = []
        startres = []

        fig=plt.figure(0)
        ax=fig.add_subplot(nrow, ncol, ib+1, adjustable='box',aspect='equal')
        ax.cla()

        try:                
            ## need to run this to get the correct likelihood.
            ##   therefore may need to adjust prior ranges

            data, xyrange = read_data_Planck_TOI([fdir+f])
                                                            
            like, prop_sigmas, start_params = setup_sampler(
                data, xyrange,
                useNormalizedBeam=useNormalizedBeam,sigminmax=sigminmax,
                rangeScale=rangeScale)
                
            if startParams is not None:
                start_params=startParams

            res[ib].append(sample1beam(like, nMC=nMC, prop_sigmas=prop_sigmas,
                                   start_params=start_params, fac=fac, rotateParams=rotateParams,
                                   noCorrelations=noCorrelations, doBlock=doBlock))

            if figName:
                fig.savefig(figf+str(fig.number).strip()+'.png')


            sys.stdout.flush()
            fig=plt.figure(1)
            ax=fig.add_subplot(nrow, ncol, ib+1)
            samples = cat([ s.samples for s in res[ib][-1][0] ])

            for var in samples.transpose(): ax.plot(var)
            if figName:
                fig.savefig(figf+str(fig.number).strip()+'.png')

            fig=plt.figure(2)
            plt.subplots_adjust(wspace=0.3, hspace=0.25)
            pidx = np.where(like.model.unpackage(prop_sigmas)>0)[0]                    
            getdist.histgrid(res[ib][-1][0][-1], params=pidx)
            plt.subplots_adjust()

            if figName:
                fig.savefig(figf+str(fig.number).strip()+'.png')

        except None:
            print "Unexpected error:", sys.exc_info()[0]
            print "... when running File: %s, iter: %d, MC: %d" % (file, it, MC)
                                

        if closeFigs: plt.close('all')

    plt.show()
    return res



def makereport(reslist, file=sys.stdout, hasRuns=False):
    """
    print out a 'report' of the MCMC run
    """

    if not hasRuns: reslist=[reslist]
    
    for irun, resrun in enumerate(reslist):
        file.write("Run: %d\n" % irun)
        for res in resrun:

            val = res[0][1][0]
            sig = res[0][1][1]
            
            for (v, s) in zip(val, sig):
                file.write("   %f +- %f" % (v, s))

            file.write("\n")

    
