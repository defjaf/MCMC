from __future__ import division
from __future__ import with_statement

import math
import string
import os.path

import pylab
import pyfits
import MCMC
import getdist
from binnedCl.binnedClLikelihood import binnedClLikelihood
from binnedCl.binnedClModel import binnedClModel, FisherWindows, plotcorrmat
from ClData import ClData
from ClData.readbins import readbins

from pylab import *

from numpy.random import uniform
from numpy import arange, array, float64, transpose, zeros, ones, exp, logical_not
from numpy import concatenate as cat
import numpy as N

filename = "data_list.txt"
#filename = "wang_dat.txt"

mapdir = 'cmb/misc-data/MAP'
homedir = os.path.expandvars('${HOME}/home')
if not os.path.exists(homedir):
    homedir = os.path.expandvars('${HOME}')
mapdir = '/'.join( (homedir, mapdir) )

def main(nMC=(1000,), gridPlot=True, testshape=True, no_pol=False, data=None, bins=None, rotateParams=True, 
         binfile=None, prefix=None):
    N.set_printoptions(precision=4, linewidth=150, suppress=True)
    mod = binnedClModel

    #mapf = '/'.join( (mapdir, 'models/cmb_04546021.fits') )
    #mapd = pyfits.getdata(mapf)
    #ll = arange(mapd.shape[0])
    #norm = 1e6
    #ClTT = array(norm**2*mapd.field(0))
    #llCl = ClTT*ll*(ll+1)/(2*math.pi)

    #max_ell = len(llCl)-1
    #max_ell = 2100

    tmp = N.fromfile("CarloClModel.dat", sep=" ")
    tmp.shape = -1, 6
    ells = int_(tmp.T[0])  ## nb doesn't usually start at l=0
    
    max_ell = max(ells)
    ll = arange(max_ell+1)
    llClTT = N.zeros(max_ell+1)
    llClEE = N.zeros(max_ell+1)
    llClTE = N.zeros(max_ell+1)
    llClTT[ells] = tmp.T[1]
    llClEE[ells] = tmp.T[2]
    llClTE[ells] = tmp.T[3]

    manybins = True
    onebin = False

    if data is None:
        data = ClData.getClData(filename, no_pol=no_pol)
    
    # print "unsetting beam and calib uncertainty"
    # for d in data:
    #  d.beam_uncertain=False
    #  d.calib_uncertainty = 0.0
    #  #d.has_xfactors=False

    if manybins:
        if binfile is None:
    #        bins = [ 2, 11, 21, 31, 41, 51, 61, 81, 101, 121, 141, 161, 181, 201,
    #                 221, 241, 261, 281, 301, 351, 401, 451, 501, 551, 601,
    #                 651, 701, 801, 901, 1001,  max_ell]
    ##        bins = [ 2, 21, 41, 61, 101, 141, 181, 221,  261, 301, 401, 501, 601,
    #        bins = [ 45, 61, 101, 141, 181, 221,  261, 301, 401, 501, 601,
    #                 701, 801, 1001, max_ell]
            binsTT = [50, 101, 151, 201, 251, 301, 351, 401, 451, 501, 551, 601, 651, 701, 
                    751, 801, 851, 901, 971, 1031, 1091, 1151, 1211, 1271, 1331, 1391, 
                    1451, 1511, 1571, 1651, 1751, 1851, 1951] #, len(llClTT)-1]
            binsTT = [101, 201, 301, 401, 501, 601, 701, 
                      801, 901, 1031, 1151, 1271, 1391, 
                      1511, 1651, 1851, 2051] #, len(llClTT)-1]
                
            #binsTT = [2,5,10,20,30,51,101,151, 251, 351, 451, 551, 651, 
            #                751, 851, 951, 1051, 1500]  
            #bins = [50, 1031, 1091, 1151, 1211, 1271, 1331,
            # bins = [50, 1391, 1451, 1511, 1571, 1651, 1751, 1851, 1951, len(llClTT)-1]
            # bins =[101, 351, 551, 651, 731, 791, 851, 911, 971, 1031, 1091, 1151, 1211, 
            #        1271, 1331, 1391, 1451, 1511, 1571, 1651, 1751, 1851, 1951, 2101, 
            #        2301, 2501, 3000]
            # bins = [50, 1851, 1951, len(llClTT)-1]
    ##        bins = [ 50, 61, 141, 221, 301, 501, 701, 1001, len(llCl)-1]
    ##        bins = [ 2, 61, 221, 501, 1001, max_ell]
            binsTE = [164, 245, 326, 407, 488, 569, 650, 731, 812, 893, 974, 
                      1055, 1136, 1217, 1298, 1379, 1460, 1540]
            binsEE = [164, 245, 326, 407, 488, 569, 650, 731, 812, 893, 974,
                      1055, 1136, 1217, 1298, 1379, 1460, 1540]
            binsTE = [164, 326, 488, 650, 812, 974, 1136, 1540]
            binsEE = [164, 326, 488, 650, 812, 974, 1136, 1540]
        
            ells = []
            if not no_pol:
                binlist = [binsTT, binsTE, binsEE]
            else:
                binlist = [binsTT]
            
            if bins is not None:
                binlist=bins

            for bins in binlist:
                for i, b in enumerate(bins[:-1]):
                    bins[i] = (b, bins[i+1]-1)  ## nb. non-pythonic: beginning and end

            bins = [b[:-1] for b in binlist]


        else:
            bins = readbins(binfile)
        
        if not no_pol:
            shapelist = array([llClTT, llClTE, llClEE])
        else:
            shapelist = array([llClTT])
        
        
        nbins = [len(b) for b in bins]
        npar = sum(nbins)

        #ell = [int((b[0]+b[1])/2) for b in bins]
        ells = [array(b).mean(axis=1) for b in bins]
        
        ## start at a reasonable model
        if not testshape:
            Clbins = [b for b in ell if b<len(llClTT) ]
            start_params = zeros(shape=(npar,), dtype=float64) + 200.0
            #start_params[0:len(Clbins)] = llClTT[Clbins]
            start_params = llClTT[int_(ell)]
            llClshape = 1.0
            prop_sigmas = zeros(npar, float64) + 100.0
        else:
            llClshape = shapelist.copy()
            start_params = ones(shape=(npar,), dtype=float64)
            prop_sigmas = zeros(shape=(npar,), dtype=float64)+0.1
            prop_sigmas[nbins[0]:]=1.0
            

    elif onebin:
        bins = [(2, len(llClTT)//2), (len(llClTT)//2+1,len(llClTT)-1)]
        npar = len(bins)
        llClshape = llClTT
        start_params = ones(shape=(npar,), dtype=float64)
        prop_sigmas = zeros(shape=(npar,), dtype=float64) + 0.5

    mod.setBinning(bins, shapefun=llClshape)

    ell = [array(b).mean(axis=1) for b in bins]

    #ell = [(b[0]+b[1])/2 for b in bins]
    print 'bins:'
    print ell
    
    print mod.nbins
    
    like = binnedClLikelihood(data=data, model=mod)
    
    pylab.figure(0)
    
    fac = 2.4/sqrt(len(ell))   ## equiv to fac=None
    #fac*=3
    
    retval = MCMC.sampler(like, nMC, prop_sigmas, start_params, plotter=plotter,
                        fac=fac, noCorrelations=True, doBlock=True, rotateParams=rotateParams)
    
    if no_pol:
        pylab.plot(ll, llClTT, hold='true')
        pylab.xlim(0,2000)
    else:
        pylab.subplot(2,2,1)
        pylab.plot(ll, llClTT, hold='true')
        pylab.xlim(0,2000)
        pylab.subplot(2,2,2)
        pylab.plot(ll, llClTE, hold='true')
        pylab.xlim(0,2000)
        pylab.subplot(2,2,3)
        pylab.plot(ll, llClEE, hold='true')
        pylab.xlim(0,2000)
        
    pylab.figure(1)
    samples = cat([ s.samples for s in retval[0] ])
    for var in samples.transpose(): pylab.plot(var)
    
    s = retval[0][-1]
    
    mean = mod.bandpowers(s.mean())
    stdv = mod.bandpowers(s.stdev())
    covar = s.covar(unNormalized=True)
    Clcovar = mod.ClCovar(covar)
    corr = s.covar()
    
    WBl = FisherWindows(Clcovar, bins=bins, isCovar=True)
    
    nTT = len(ell[0])
    nTE = len(ell[1])
    nEE = len(ell[2])
    if prefix is not None:
        for es, ms, ss, suf in zip(ell, mean, stdv, [".bp", ".bpte", ".bpee"]):            
            with open(prefix+suf, "w") as f:
                for ell1, mean1, stdv1 in zip(es, ms, ss):
                    print >> f, ell1, mean1, stdv1, stdv1
                
        with open(prefix+".covar", "w") as f:
            N.savetxt(f, Clcovar, fmt="%f")
            
        with open(prefix+".corr", "w") as f:
            N.savetxt(f, corr, fmt="%f")
            
        for ibin, win in enumerate(WBl):
            with open(prefix+str(ibin+1), "w") as f:
                for l, Wl in enumerate(win.T):
                    print >> f, l, Wl[0], Wl[1], Wl[2]
                
    if gridPlot:
        pylab.figure(2)
        getdist.histgrid(s)
    else:
        getdist.printvals(s)
    
    # params=xrange(s.samples.shape[1])  
    # mean = N.empty(s.samples.shape[1])
    # stdv = N.empty(s.samples.shape[1])
    # for param in params:
    #     s1 = s.samples.T[param]
    #     mean[param] = s1.mean()
    #     stdv[param] = s1.std()
    
    

    for l, m, e in zip(ell, mean, stdv):
        for l1, m1, s1 in zip(l, m, e):
            print '%d %f %f' % (int(l1),m1,s1)       

    ## split this off into mod.plotWin?
    nspec = len(bins)
    iwin = 0
    for ifig, binlist in enumerate(bins):  ## loop over figures
        pylab.figure(5+ifig)
        nr = nc = int(sqrt(len(binlist)))
        if nr*nc<len(binlist): nr+=1
        if nr*nc<len(binlist): nc+=1
        for ipanel, bin in enumerate(binlist): ## loop over panels
            pylab.subplot(nr, nc, ipanel+1)
            for ispec in xrange(nspec):  ## loop over spectrum types 
                plot(WBl[iwin, ispec])
            axvspan(bin[0], bin[1],  facecolor='g', alpha=0.25)
            axhline(0, color='k')
            iwin += 1
            
    plotcorrmat(corr, bins=[nTT, nTE, nEE])
    
    return retval


def plotter(sampler):
    """
    make a plot from the results of the sampler
    (could probably make this more generic if plotMod is a member fn of the model)
    for now, only plot the *first* dataset (since for this case that's the 'good' data)
    """
    
    mod = sampler.like.model
    data = sampler.like.data

    params = where(sampler.prop.sigmas>0)[0]


    ntot = len(sampler.samples)
    stride = 1 #max(1,ntot//sampler.naccept)
    ana = MCMC.chain_analyze(sampler.samples[(ntot//5)::stride,:], params=params)
    vals = sampler.like.model.package(ana[0])
    sigs = sampler.like.model.package(ana[1])
    
    ### need to take shape into account
    vals = mod.bandpowers(vals)
    sigs = mod.bandpowers(sigs)
    
    
    #vals = vals*mod.ellctr*(mod.ellctr+1)/(2*math.pi)
    #sigs = sigs*mod.ellctr*(mod.ellctr+1)/(2*math.pi)


    print vals  #sampler.like.model.fmtstring % tuple(ana[0])
    print sigs  #sampler.like.model.fmtstring % tuple(ana[1])
    print ana[2]
    

    ### or replace with mod.plotmod if written...
    ## nb. bandpowers(vals) makes into 3*n list if polarized
    nrc = 1
    if len(vals)>1: nrc=2
    if mod.nbins != mod.nparam:
        for iplot, (l, v, s) in enumerate(zip(mod.ellctr, vals, sigs)):
            pylab.subplot(nrc,nrc,iplot+1)
            pylab.cla()
            pylab.errorbar(l, v, yerr=s, fmt='bo')
    else:
        pylab.cla()
        pylab.errorbar(mod.ellctr, vals, yerr=sigs, fmt='bo')
        
    #m = mod(vals); c1 = m()[0]; ell = N.arange(c1.size)
    #c = c1*ell*(ell+1)/(2*N.pi)
    #pylab.plot(ell, c)
    

########################################################
def getlike(ibins=[1], data=None):
    mod = binnedClModel

    #mapf = '/'.join( (mapdir, 'models/cmb_04546021.fits') )
    #mapd = pyfits.getdata(mapf)
    #ll = arange(mapd.shape[0])
    #norm = 1e6
    #ClTT = array(norm**2*mapd.field(0))
    #llCl = ClTT*ll*(ll+1)/(2*math.pi)

    #max_ell = len(llCl)-1
    #max_ell = 2100

    tmp = N.fromfile("CarloClModel.dat", sep=" ")
    tmp.shape = -1, 6
    ells = int_(tmp.T[0])  ## nb doesn't usually start at l=0

    max_ell = max(ells)
    ll = arange(max_ell+1)
    llClTT = N.zeros(max_ell+1)
    llClEE = N.zeros(max_ell+1)
    llClTE = N.zeros(max_ell+1)
    llClTT[ells] = tmp.T[1]
    llClEE[ells] = tmp.T[2]
    llClTE[ells] = tmp.T[3]

    manybins = True
    onebin = False

    testshape = True

    if data is None: data = ClData.getClData(filename, no_pol=True)

    # print "unsetting beam and calib uncertainty"
    # for d in data:
    #     d.beam_uncertain=False
    #     d.calib_uncertainty = 0.0
    #     #d.has_xfactors=False

    if manybins:
#        bins = [ 2, 11, 21, 31, 41, 51, 61, 81, 101, 121, 141, 161, 181, 201,
#                 221, 241, 261, 281, 301, 351, 401, 451, 501, 551, 601,
#                 651, 701, 801, 901, 1001,  max_ell]
##        bins = [ 2, 21, 41, 61, 101, 141, 181, 221,  261, 301, 401, 501, 601,
#        bins = [ 45, 61, 101, 141, 181, 221,  261, 301, 401, 501, 601,
#                 701, 801, 1001, max_ell]
        bins = [50, 101, 151, 201, 251, 301, 351, 401, 451, 501, 551, 601, 651, 701, 
                751, 801, 851, 901, 971, 1031, 1091, 1151, 1211, 1271, 1331, 1391, 
                1451, 1511, 1571, 1651, 1751, 1851, 1950, len(llClTT)-1]
        bins = [50, 1651, 1751, 1851, 1951, 3000] #, len(llClTT)-1]
        bins = [1, 1851, 1951, len(llClTT)-1]
        
##        bins = [ 50, 61, 141, 221, 301, 501, 701, 1001, len(llCl)-1]
##        bins = [ 2, 61, 221, 501, 1001, max_ell]
        for i, b in enumerate(bins[:-1]):
            bins[i] = (b, bins[i+1]-1)  ## nb. non-pythonic: beginning and end
        bins = bins[:-1]
        npar = len(bins)

        ell = [int((b[0]+b[1])/2) for b in bins]
        ## start at a reasonable model
        if not testshape:
            Clbins = [b for b in ell if b<len(llClTT) ]
            #start_params = zeros(shape=(npar,), dtype=float64)
            start_params = llClTT[int_(ell)]
            shape = 100.
            prop_sigmas = zeros(npar, float64) + 100.0
        else:
            shape = llClTT
            start_params = ones(shape=(npar,), dtype=float64)
            prop_sigmas = zeros(shape=(npar,), dtype=float64) + 0.1


    elif onebin:
        bins = [(2, len(llClTT)//2), (len(llClTT)//2+1,len(llClTT)-1)]
        npar = len(bins)
        shape = llClTT
        start_params = ones(shape=(npar,), dtype=float64)
        prop_sigmas = zeros(shape=(npar,), dtype=float64) + 0.5
    
    print 'Start:'
    print start_params

    mod.setBinning(bins, shapefun=shape)

    ell = [(b[0]+b[1])/2 for b in bins]
    print 'bins:'
    print ell
    
    like = binnedClLikelihood(data=data, model=mod)
    
    #### plot the likelihood as a function of power in a single bin
    for ibin in ibins:
        bps = N.linspace(0, 5*start_params[ibin], 100)
        lnlikearr = N.zeros_like(bps)
        for i, bp in enumerate(bps):
            pars = start_params.copy()
            pars[ibin] = bp
            lnlikearr[i] =  like.lnLike(pars)
            
        lnlikearr[logical_not(isfinite(lnlikearr))] = min(lnlikearr[isfinite(lnlikearr)])-10.0
        lnlikearr = lnlikearr-max(lnlikearr)
        plot(bps, lnlikearr)

        likearr = exp(lnlikearr)

        avg = (bps*likearr).sum()/likearr.sum()
        print '<par>=%f' % (avg)
        
    
    return data, like, likearr
    
    
########################################################

import hotshot, hotshot.stats
def profrun():
    #profile.Profile.bias=1e-5
    #profile.run('test_binnedCl.main(nMC=100)', 'Clprof.out')
    prof = hotshot.Profile("binnedCl.prof")
    profout = prof.runcall(main, nMC=100)
    prof.close()
    return profout
