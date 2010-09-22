#!/usr/bin/env python
# encoding: utf-8
"""
driver.py

Created by Andrew H. Jaffe on 2010-05-22.
Copyright (c) 2010 Imperial College London. All rights reserved.
"""

from __future__ import division

import sys
import os

import operator
import cPickle as pickle

import numpy as np

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("AGG")
    
import matplotlib.pyplot as plt

import MCMC
import likelihood
import data
import model
import getdist

import joblib


#### pattern after test_binnedCl.py and/or BeamFit/driver.py BeamFit/MAXIPOLBeamData.py

fname_DLC = "./submmSED.txt"
fname_ERCSC = "./submmSED/ercsc_iifscz.txt"
fname_MRR = "./submmSED/ERCSCalliifscz4550850.dat"


### TODO: add calculation of the likelihood/posterior of the posterior mean params

def main(filename=fname_MRR, i=0, rotateParams=False, onecomponent=True, getNorm=True, start=None, sigmas=None, 
         nMC=(10000,10000), nDerived=None, noPlots=False, DLC=False, MRR=True, fig0=0, savefig=False, retMCMC=True,
         fdir = "./"):
        
        
    ret = []
        
    if not operator.isSequenceType(i):
        idata = [i]
    else:
        idata = i
        
    ## read the data
    if DLC:
        alldata = data.readfluxes_DLC(filename)
    elif MRR:
        alldata = data.readfluxes_MRR(filename)
    else:
        alldata = data.readfluxes_ERCSC_TopCat(filename)
        
    for i in idata:
        
        try:
            dat = alldata[i]
            name = dat.name + "(z=%4.2f)" % (dat.z)
        except TypeError:
            dat = [alldata[j] for j in i]
            name = " + ".join(str(int(d.name)) for d in dat)
        
        print "Object[s] %s" % name
    
        ## initialize the model (with a class, not an object)
        if getNorm:
            if onecomponent:
                mod = model.submmModel1_normalized_logA
            else:
                mod = model.submmModel2_normalized_logA
            like = likelihood.SEDLikelihood_normalized(data=dat, model=mod)        
        else:
            if onecomponent:
                mod = model.submmModel1
                like = likelihood.SEDLikelihood1(data=dat, model=mod)
            else:
                mod = model.submmModel2            
                like = likelihood.SEDLikelihood2(data=dat, model=mod)
        
        if start is None: 
            start_params = np.asarray(mod.startfrom(random=False))
        else:
            start_params = np.asarray(start)
        if sigmas is None:
            prop_sigmas = start_params/4.0
        else:
            prop_sigmas = sigmas
        
        if nDerived is not None:
            like.nDerived = nDerived
    
        mcmc, ana = MCMC.sampler(like, nMC, prop_sigmas, start_params, plotter=None,
                            fac=None, noCorrelations=True, doBlock=True, rotateParams=rotateParams)

        if not noPlots:

            fig = plt.figure(fig0)
            lnLike = []
            maxlnLike, maxLikeParams = getdist.histgrid(mcmc[-1])
            plt.suptitle(name)
    
            if savefig:
                try:
                    fname = fdir+name + savefig
                except TypeError:
                    fname = fdir+name
                    
                fig.savefig(fname+"_0.png")
                plt.close(fig)
    
            fig=plt.figure(fig0+1)
            params = ana[0]
            meanmod = mod(*params)
            meanlnProb = like.lnLike(params)
            print "ln Pr of mean = %f" % meanlnProb
            MLmod = mod(*maxLikeParams)
            try:
                meanmod.plot(dat, wavelength=True, logplot=True)
                MLmod.plot(dat, wavelength=True, logplot=True)
                plt.suptitle(name)    
            except AttributeError:
                pass
            if savefig:
                fig.savefig(fname+"_1.png")
                plt.close(fig)
                
            fig0 += 2
    
        if retMCMC:
            ret.append((mcmc, ana, (maxlnLike, maxLikeParams, meanlnProb), name))  ### meanlnProb added
        else:
            del mcmc
            ret.append((ana, (maxlnLike, maxLikeParams, meanlnProb), name))
    
    if len(ret) == 1:
        ret = ret[0]
        
    return ret
    
    

# idata =[i*50 for i in range(6,29)]
# nMC = (15000,100000)

#idata =[i*25 for i in range(12,57)]
nMC = (15000,100000)
#    fil = "./ercsc_iifscz.txt"
idata = range(0,726)
fil = "./ERCSCalliifscz4550850.dat"


def many(which = range(3), idata=idata, nMC = nMC, fil=fil, fdir="./"):

    print "Using file %s" % fil
    
    ret1 = ret2 = ret3 = []

    if 0 in which:
        print "Two-Component beta = 2"
        ret1 = main(fil, getNorm=True, i = idata, 
                    start=(1,2.,10,0.1,2.,20), sigmas=(1,0,2, 1, 0, 2), retMCMC=False,
                    nMC=nMC, onecomponent=False, fig0=0, savefig="_2comp_b2", fdir=fdir)

    if 1 in which:
        print "One-Component"
        ret2 = main(fil, getNorm=True, i = idata, 
                    start=(1,2.,10), sigmas=(1,2,2), retMCMC=False,
                    nMC=nMC, onecomponent=True, fig0=100, savefig="_1comp", fdir=fdir)
                
    if 2 in which:
        print "One-Component beta = 2"
        ret3 = main(fil, getNorm=True, i = idata, 
                    start=(1,2.,10), sigmas=(1,0,2), retMCMC=False,
                    nMC=nMC, onecomponent=True, fig0=200, savefig="_1comp_b2", fdir=fdir)
                
    return ret1, ret2, ret3
    

def postprocess(dirname="./"):
    
    
    ret = []
    idxs = [[0,2,3,5], [0,1,2], [0,2]]
    for i in range(3):
        
        fname = dirname+"out_[%d].pickle" % i        
        try:
            with open(fname) as f:
                ret0 = pickle.load(f)
        except IOError:
            continue
            
        ix = idxs[i]
        ret0 = ret0[i]
        nobj = len(ret0)
        npar = len(ix)
        
        dt = np.dtype([
            ('name', 'S19'),
            ('mean', np.float, (npar,)), 
            ('sig', np.float, (npar,)), 
            ('covar', np.float, (npar,npar)), 
            ('ML', np.float),
            ('ev', np.float),
            ('MLpar', np.float, (npar,)),
            ('MeanL', np.float),
            ('evMean', np.float),
            ('dlnLike', np.float)
        ])
        ret_i = np.empty(nobj, dtype = dt)
        
        ### each of ret[i] has format
        ## [ <for each object>... 
        ##   (
        ##     ([parameter means], [parameter sigmas], [normalized covar]),
        ##     (ML, [ML params], meanL),  !!! meanL added recently
        ##     name    !!! added in recent versions
        #    )
        ## ...]
 
        for iobj, obj in enumerate(ret0):
            ret_i[iobj]['name'] = obj[2]
            ret_i[iobj]['mean'][:] = np.array(obj[0][0])[ix]
            ret_i[iobj]['sig'][:] = np.array(obj[0][1])[ix]
            covar = np.array(obj[0][2])[:,ix][ix]
            ML = obj[1][0]
            ret_i[iobj]['covar'][:,:] = covar
            ret_i[iobj]['ML'] = ML
            ret_i[iobj]['ev'] = ML + 0.5*np.linalg.det(covar) + npar*0.5*np.log(2*np.pi)
            ret_i[iobj]['MLpar'][:] = obj[1][1][ix]
            try:
                meanL = obj[1][2]
                ret_i[iobj]['ML'] = meanL
                ret_i[iobj]['evMean'] = meanL + 0.5*np.linalg.det(covar) + npar*0.5*np.log(2*np.pi)
                ret_i[iobj]['dlnLike'] = ML-meanL
            except IndexError:
                pass
        
        ret.append(ret_i)
    
    return ret


def writeTabAll(ret123, fbase, ext='.npy', dirname=None):
    if dirname is not None:
        ret123 = postprocess(dirname)
        
    for i, r in enumerate(ret123):
        fname = fbase + str(i) + ext
        writeTab(r, fname)
        
    return ret123
        

def writeTab(ret, fname, names=None):
    """ write the output of the postprocess function to a text file 
        run separately on each of the elements of postprocess()
     """
     
    try:
        anames = ret['name']
    except ValueError:
        anames = np.array(names)

    alls = np.hstack([anames.reshape(-1,1), ret['MLpar'], ret['sig'], ret['dlnLike'].reshape(-1,1)])
        
    npar = len(ret['MLpar'][0])
    hdr = ['Name']
    for i in range(npar):
        hdr.append("ML param %d" % i)
    for i in range(npar):
        hdr.append("sigma param %d" % i)
    hdr.append("dlnLike")
        
    hdr = ("%20s "*(2*npar+2)) % tuple(hdr)
    with open(fname, 'w') as f:
        f.write(hdr + '\n')
        np.savetxt(f, alls, fmt="%20s", delimiter=' ')
        
        
def plotter(sampler):
    """
    make a plot from the results of the sampler
    (could probably make this more generic if plotMod is a member fn of the model)
    for now, only plot the *first* dataset (since for this case that's the 'good' data)
    
    TODO: calculate stats for the amplitude here?
    
    NOT USED IN ABOVE YET
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


    #m = mod(vals); c1 = m()[0]; ell = N.arange(c1.size)
    #c = c1*ell*(ell+1)/(2*N.pi)
    #pylab.plot(ell, c)
    
def simul():
    """ supposed to do embarassingly parallel python launches, but doesn't appear to work... """
    par = joblib.Parallel(n_jobs=3, verbose=1)
    ret123 = par(joblib.delayed(many([w])) for w in range(3))
    with open("out_123.pickle", 'w') as f:
        pickle.dump(ret123, f)
        
    return ret123


if __name__ == '__main__':
    fdir = "./figs0921a/"
    which = []
    for s in reversed(sys.argv):
        try:
            which.append(int(s))
        except ValueError:
            break
    print "which=", which
    ret = many(which, fdir=fdir)
    with open("out_"+"".join(str(which).split(' '))+".pickle", 'w') as f:
        pickle.dump(ret, f)
    
    

