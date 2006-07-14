from __future__ import division

import math
import string
import os.path
import cPickle

import pylab
import MCMC
import getdist

import topo_model
from topo_likelihood import topo_likelihood
from likelihood.likico import likico
from likelihood.likoct import likoct
from likelihood.likdih import likdih
from likelihood.liktetr import liktetr

from numpy.random import uniform
from numpy import arange, array, float64, transpose, zeros, ones, where, nonzero
from numpy import concatenate as cat
import numpy as N

    
def like_grid(modname='ico', almfile=None):
    """
    compute n-dimensional grid of topology likelihood
    (could make much more generic?)
    
    """
    N.set_printoptions(precision=4, linewidth=150, suppress=True)
    
    lik = { 'ico': likico, 'oct': likoct, 'dih': likdih, 'tetr': liktetr}

    datdir = "./topology/likelihood/lik"+modname+"/dat/\0"
    if almfile is None:
        # almfile = "alm64_1.dat\0"
        #    almfile = "wmapalm.dat\0"
        almfile = "alm3yrall.dat\0"

    print 'using datdir=%s' % datdir
    likfun = lik[modname].alikelihood
    
    lik[modname].readdata(N.array(datdir), N.array(almfile))

    ang_lims = {'ico':  [2/5, 1, 2/5],
                'oct':  [1/2, 1, 1/2],
                'dih':  [1/2, 1, 1/2], 
                'tetr': [1/2, 1, 1/2]}

    for key, lims in ang_lims.iteritems():
        ang_lims[key] = N.array(lims) * math.pi

    al = ang_lims[modname]
### amplitude, alpha, beta, gamma, H0        
    param_min = (5.0e-5, 54.0, 0.0,   0.0,   0.1)
    param_max = (5.0e-3, 72.0, al[0], al[1], al[2])
    nstep =     (5,)*5
    
    npar = len(nstep)
    assert len(param_min)==len(param_max)==npar, 'Bad number of parameters'

    #kind of inelegant, but it does the job...
    args = tuple(slice(p1, p2, n*1j) for p1,p2,n in zip(param_min, param_max, nstep))
    param_steps = N.ogrid.__getitem__(args)
    param_grid = N.mgrid.__getitem__(args)
    param_list = N.transpose(N.asarray(param_grid).reshape(len(param_min), -1))

    #for p1, p2, n in zip(param_min, param_max, nstep):
    #    params.append(N.linspace(p1, p2, n))

    like = N.empty(nstep, dtype=float64)

    ## now need to iterate over all permutations of the params grid

    viewstep = 100
    for ip, par in enumerate(param_list):
        if not (ip % viewstep): print ip, par, ' ',
        try:
            like1 = likfun(*par)
        except:
            print "\nlikelihood failure at ", ip, par
            like1 = N.nan
            
        like.flat[ip] = like1
        if not (ip % viewstep): print like1

    return param_grid, like


def analyze_grid(like):
    """
    analyze an n-dimensional likelihood grid
    (calculate means, marginalized likelihoods, find maxima)
    """

    
#### want to make an iterator which takes a sequence of sequences (!)
####    and returns a tuple 
