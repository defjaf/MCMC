from __future__ import division

import math
import string
import os.path
import cPickle

import pylab
import MCMC
import getdist_ahj as getdist

import topo_model
from topo_likelihood import topo_likelihood
from likelihood.likico import likico
from likelihood.likoct import likoct
from likelihood.likdih import likdih
from likelihood.liktetr import liktetr

from numpy.random import uniform
from numpy import arange, array, float64, transpose, zeros, ones, where, nonzero
from numpy import concatenate as cat
import numpy

noplot = False

def testall():
    test('ico')
    test('oct')
    test('dih')
    test('tetr')


def test(modname='ico', almfile=None):
    
    lik = { 'ico': likico, 'oct': likoct, 'dih': likdih, 'tetr': liktetr}
    
    datdir = "./topology/likelihood/lik"+modname+"/dat/\0"
    if almfile is None:
        # almfile = "alm64_1.dat\0"
        #    almfile = "wmapalm.dat\0"
        almfile = "alm3yrall.dat\0"
    
    print 'using datdir=%s' % datdir
    likfun = lik[modname].alikelihood
    
    lik[modname].readdata(numpy.array(datdir), numpy.array(almfile))
    
    alpha=0. #0.364298
    beta=0. #0.1290
    gamma=0.# 1.84952
    hmin = 52.
    hmax = 70.
    hub = 64.0
    dh = 1.0
    nstep = int(1+(hmax - hmin)/dh)
    a = []; h=[]; l=[]
    for i, hub in enumerate(numpy.arange(hmin, hmax, dh)):
        #hub=dh*(hmax-hmin)/(nstep-1)*i+hmin
        #alpha = math.pi/nstep*i
        try:
            like = likfun(1.0,hub,alpha,beta,gamma)
        except:
            print 'likelihood failure at alpha, hub=', alpha, hub
        a+= [alpha]; h+=[hub]; l+=[like]
        print  alpha, hub, like
    
    if not noplot:
        pylab.plot(h,l)
    
    return a, h, l
    


def main(nMC=(100, 300, 1000), noCorrelations=True, fac=None, doBlock=True,
         almfile=None, topo='ico', fig=0, doSim=False):
    
    numpy.set_printoptions(precision=4, linewidth=150, suppress=True)
    
    datdir = "./topology/likelihood/lik"+topo.strip()+"/dat/\0"
    if almfile is None:
    #   almfile = "alm60_1.dat\0"
    #   almfile = "alm64_1.dat\0"
    #    almfile = "wmapalm.dat\0"
       almfile = 'alm3yrall.dat'
    mod = topo_model.topo_model
    like = topo_likelihood(model=mod, datdir=datdir, almfile=almfile,
                           topo=topo)
    npar = 5
    
    print "Using likelihood: %s" % like.topo
    
    doBlock=True
    
    mod.H0_min = 54.0
    mod.H0_max = 68.0
    
    H0_start = uniform(mod.H0_min, mod.H0_max)   #64
    
    start_params = (1.0e-3, (0.0,0.0,0.0), H0_start)
    #prop_sigmas = (0.3e-3, (0.05, 0.05, 0.05), 3.0)
#    prop_sigmas = (3.0e-3, (0.5, 0.5, 0.5), 5.0)
    prop_sigmas = (3.0e-3, (2.5, 2.5, 2.5), 5.0)
    mod.paramBlocks = [0,1,2,3,4]
    #mod.paramBlocks = [-1, 0,1,2,3]
    mod.nBlock = 5
    
    if doSim:
        print "Assuming simulation amplitude ~ 1"
        start_params=(1.0, (0.0,0.0,0.0), H0_start)
        prop_sigmas = (0.3, (0.5, 0.5, 0.5), 2.0)
    
    print ("Starting point:  " + mod.fmtstring) % tuple(mod.unpackage(start_params))
    print ("Starting sigmas: " + mod.fmtstring) % tuple(mod.unpackage(prop_sigmas))
    
    res, ana = MCMC.sampler(like, nMC, prop_sigmas, start_params,
                        plotter=plotter, fac=fac,
                        noCorrelations=noCorrelations, doBlock=doBlock)
    
    samples = cat([ s.samples for s in res ])
    if not noplot:
        pylab.figure(10*fig)
        for var in (samples.transpose())[1:]: pylab.plot(var)
        pylab.show()
        
        pylab.figure(10*fig+1)
        getdist.histgrid(res[-1], [0,1,2,3,4])
        pylab.show()
    
    return res, ana

def plotter(sampler):
    """
    show results; don't plot...
    """
    
    mod = sampler.like.model
    data = sampler.like.data
    
    params = nonzero(sampler.prop.sigmas>0)[0]  ### [0] added since tuple returned
    
    ntot = len(sampler.samples)
    stride = 1 #max(1,ntot//sampler.naccept)
    ana = MCMC.chain_analyze(sampler.samples[(ntot//5)::stride,:], params=params)
    vals = sampler.like.model.package(ana[0])
    sigs = sampler.like.model.package(ana[1])
    print vals  #sampler.like.model.fmtstring % tuple(ana[0])
    print sigs  #sampler.like.model.fmtstring % tuple(ana[1])
    print ana[2]

def doall(nMC=(100,500,10000), file='wmap06res_10000.pickle'):
    
    res = {}
    for i, x in enumerate(['ico', 'oct', 'tetr', 'dih']):
        res[x] = main(almfile='alm3yrall.dat\0', nMC=nMC, topo=x, fig=10*i)
    
    pkl = open(file, 'w')
    cPickle.dump(res, pkl);
    pkl.close()
    
    return res

def plotall(picklefile, pre='', suff='', save=True):
    
    fp = open(picklefile)
    allres = cPickle.load(fp)
    fp.close()
    
    for i, (name, res)  in enumerate(allres.iteritems()):
        pylab.figure(i)
        print name
        getdist.histgrid(res[0][-1])
        if save:
            pylab.savefig(pre+name+suff+'.png')

    
    
