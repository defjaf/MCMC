from __future__ import division

import math

from likelihood.likico import likico
from likelihood.likoct import likoct
from likelihood.likdih import likdih
from likelihood.liktetr import liktetr

import numpy as N
from numpy import float64

def do_all():
    
    ret = {}
    dir = 'topology/out/'
    for mod in ['ico', 'oct', 'dih', 'tetr']:
        print '*********** running model %s *****************' % mod
        ret[mod] = like_grid(modname=mod, outfile=dir+'lik_'+mod+'_grid.out', nstep=12)
    
    return ret


def like_grid(modname='ico', almfile=None, datname=None, outfile=None, nstep=10):
    """
    compute n-dimensional grid of topology likelihood
    (could make much more generic?)
    
    """
    N.set_printoptions(precision=4, linewidth=150, suppress=True)
    
    if outfile is not None:
        try:
            fp = open(outfile, 'wa')
            needs_close=True
        except:
            ### assume outfile is a file-like object already; should really check!
            fp = outfile
            needs_close=False
    
    lik = { 'ico': likico, 'oct': likoct, 'dih': likdih, 'tetr': liktetr}
    
    if datname is None: datname = modname
    
    datdir = "./topology/likelihood/lik"+datname+"/dat/\0"
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
    param_min = (5.0e-5, 54.0, 0.0,   0.0,   0.0)
    param_max = (5.0e-3, 72.0, al[0], al[1], al[2])
    npar = len(param_min)
    
    if N.size(nstep) == 1: nstep = (nstep,)*npar
    
    parfmt = '%f '*npar
    
    assert len(nstep)==len(param_max)==npar, 'Bad number of parameters'
    
    #kind of inelegant, but it does the job...
    args = tuple(slice(p1, p2, n*1j) for p1,p2,n in zip(param_min, param_max, nstep))
    param_steps = N.ogrid.__getitem__(args)
    param_grid = N.mgrid.__getitem__(args)
    param_list = N.transpose(N.asarray(param_grid).reshape(npar, -1))
    
    #for p1, p2, n in zip(param_min, param_max, nstep):
    #    params.append(N.linspace(p1, p2, n))
    
    like = N.empty(nstep, dtype=float64)
    
    ## now need to iterate over all permutations of the params grid
    
    viewstep = 100
    for ip, par in enumerate(param_list):
        if not (ip % viewstep): print ip, par, ' ',
        like1 = likfun(*par)
        ## the following is a problem since it catches things like ctrl-C!
        
        #try:
        #    like1 = likfun(*par)
        #except:
        #    print "\nlikelihood failure at ", ip, par
        #    like1 = N.nan
        
        like.flat[ip] = like1
        if not (ip % viewstep): print like1
        
        if outfile is not None:
            fp.write(('%d '+parfmt+'%f\n') % ((ip,) + tuple(par) + (like1,)))
    
    if needs_close: fp.close()
    
    return param_grid, like


def analyze_grid(param_grid, lnlike, names=None, printstats=True):
    """
    analyze an n-dimensional likelihood grid
    (calculate means, marginalized likelihoods, find maxima)
    
    TODO: allow an 'ogrid' of parameters, instead of an 'mgrid'
    """
    
    assert param_grid[0].shape==like.shape, 'Bad shapes'
    
    npar = param_grid.shape[0]
    
    ## get the steps and 'deltas' in each dimension from the mgrid param_grid
    ## (i.e. recreate the ogrid)
    deltas = N.empty(npar, dtype=float64)
    steps = []
    for i in xrange(npar):
        steps.append(N.array([param_grid[(i,) + (0,)*i + (j,) + (0,)*(npar-1-i)]
                            for j in xrange(param_grid.shape[i+1])]))
        deltas[i] = steps[-1][1] - steps[-1][0]
    
    dNpar = deltas.prod()
    
    means = N.empty(npar, dtype=float64)
    covar = N.empty((npar, npar), dtype=float64)
    like = N.exp(lnlike)
    weights = like.flat
    
    lnEv = N.log(weights.sum()) + N.log(dNpar)
    
    for i in xrange(npar):
        means[i] = N.average(param_grid[i].flat, weights=weights)
    
    for i in xrange(npar):
        for j in xrange(npar):
            covar[i,j] = (N.average((param_grid[i]*param_grid[j]).flat,
                                    weights=weights)
                          -means[i]*means[j])
    
    stdevs = N.sqrt(covar.diagonal())
    
    ## a list since there aren't necessarilly the same number in each dimension
    ## could use an object array?
    ## need to iterate over dimensions
    ## assumes linear parameter spacing!
    lnlike1d = []
    for ipar, istep in enumerate(steps):
        lnlike1d.append
    
    if printstats:
        print 'ln Evidence: ', lnEv
        print means
        print stdevs
        print covar
    
    return lnEv, means, stdevs, covar
