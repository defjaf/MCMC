#!/usr/bin/env python
# encoding: utf-8
"""
M31driver.py

Created by Andrew H. Jaffe on 2012-08-26.

"""

from __future__ import division

import sys
import os

import getopt

import operator
import cPickle as pickle

import numpy as np

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("AGG")
    plottype = "png"
    
plottype = "pdf"

    #print "NOT USING LaTeX"; matplotlib.rc('text', usetex=False)  ## memory leak with latex????
    
import matplotlib.pyplot as plt

import MCMC
import likelihood
import data
import model
import M31model
from .. import getdist_ahj ##as getdist


### TODO: add calculation of the likelihood/posterior of the posterior mean params

## added random and randomrestart params

def M31(rotateParams=False, start=None, sigmas=None, 
         nMC=(10000,10000), noPlots=False, fig0=0, savefig=False, retMCMC=True,
         fdir = "./", logplot=True, check=None, wavelength=False, doBlock=True,
         filename="M31/M31Flux-DX9d.dat",
         random=False, randomrestart=False):
        
        
    speedOfLight = 299792.  ## micron GHz
    
    noHist = False  
    
    ret = []
        
#    start = np.array(start) + 0.5*np.array(sigmas)*np.random.randn(len(start))
#    print "randomizing start:", start
            
    alldata = data.readfluxes_M31(filename)
        
    dat = alldata[0]
    name = dat.name 
        
    print "Object %s" % name
    print "N_data = %d" % dat.n
    
    mod = M31model.M31model
    like = likelihood.SEDLikelihood_normalized(data=dat, model=mod)
    start_params = np.asarray(mod.startfrom(random=random))
    if sigmas is None:
        sigmas = prop_sigmas = np.abs(start_params)/4.0
    else:
        prop_sigmas = sigmas
    
    mcmc, ana = MCMC.sampler(like, nMC, prop_sigmas, start_params, plotter=None,
                        fac=None, noCorrelations=True, doBlock=doBlock, rotateParams=rotateParams,
                        randomrestart=randomrestart)

    if not noPlots:

        fig = plt.figure(fig0)
        lnLike = []
        plt.subplots_adjust(wspace=0.3, hspace=0.25)
        maxlnLike, maxLikeParams = getdist_ahj.histgrid(
                                     mcmc[-1], noPlot=noHist, params=np.where(np.array(sigmas)>0)[0],
                                     burn=0.2, stride=1)
        plt.subplots_adjust()

        #plt.suptitle(name)

        if savefig:
            try:
                fname = fdir+name + savefig
            except TypeError:
                fname = fdir+name
                
            try:     ### move close to finally clause? 
                fig.savefig(fname+"_0."+plottype)
                plt.close(fig)
            except Exception as err:  ## dangerous -- catch anything!
                print "CAUGHT Exception!! -- WHILE SAVING %s" % fname+"_0.png"
                print "Error Info: "
                print type(err)     # the exception instance
                print err.args      # arguments stored in .args
                print err

        fig=plt.figure(fig0+1)
        params = ana[0]
        meanmod = mod(*params)
        meanlnProb = like.lnLike(params)
        mean_chi2 = like.chi2(params)
        print "ln Pr of mean = %f" % meanlnProb
        print "chi2(mean) = %f"% mean_chi2
        
        MLmod = mod(*maxLikeParams)
        try:
            meanmod.plot(dat, wavelength=wavelength, logplot=logplot, components=False)
            MLmod.plot(dat, wavelength=wavelength, logplot=logplot)
            # plt.suptitle(name)    
        except AttributeError:
            pass
        if savefig:
            try:     ### move close to finally clause? 
                fig.savefig(fname+"_1."+plottype)
                plt.close(fig)
            except Exception as err:   ## dangerous -- catch anything!
                print "CAUGHT Exception!! -- WHILE SAVING %s" % fname+"_1.png"
                print "Error Info: "
                print type(err)     # the exception instance
                print err.args      # arguments stored in .args
                print err

        fig0 += 2

    ret= (ana, (maxlnLike, maxLikeParams, meanlnProb), name) ### meanlnProb added

    ### collect further information to return
    ret += (zip(dat.d, dat.sig),)
    #ret += (MLmod.flux(nu1, nu2),)
    ret += (0,)


    if retMCMC:
        ret = (mcmc,)+ret
    else:
        ret = ret
                
    if check is not None:
        with open(check, 'w') as f:
            pickle.dump(ret, f)
            
    if len(ret) == 1:
        ret = ret[0]
        
    return ret
    
    