

import math
from math import log

import numpy as N

from . import Likelihood
from .likelihood.likico import likico
from .likelihood.likoct import likoct
from .likelihood.likdih import likdih
from .likelihood.liktetr import liktetr

likfuns = {'ico': likico, 'oct': likoct, 'dih': likdih, 'tetr': liktetr}

### alpha, beta, gamma limits
ang_lims = {'ico':  [2/5, 1, 2/5],
            'oct':  [1/2, 1, 1/2],
            'dih':  [1/2, 1, 1/2], 
            'tetr': [1/2, 1, 1/2]}

for key, lims in ang_lims.items():
    ang_lims[key] = N.array(lims) * math.pi

class topo_likelihood(Likelihood.Likelihood):
    """
    represents the likelihood function for a topologically nontrivial universe
    P(alm | A, H0, alpha, beta, gamma) where
       A is propto the mean-square amplitude
       H0 controls the *curvature* of the Universe
       (alpha, beta, gamma) are the euler angles controlling the orientation of the
            fundamental domain
    """

    def __init__(self, data=None, model=None, datdir="", almfile="",
                 topo='ico'):
        self.set_topo(topo)
        print("reading likelihood files (data=%s)... "%almfile, end=' ')
        
        self.likfun.readdata(datdir, almfile)
        print("done")
        Likelihood.Likelihood.__init__(self,data, model)
        self.model.ang_lims = ang_lims[topo]
        

    def lnNorm1(self):
        return 0.0
    
    def setModel(self, params):
        """set the model for the given parameters"""
        self.this_model = self.model(*params)
        self.model_vals = self.this_model()

    def set_topo(self, topo):
        self.likfun = likfuns[topo]
        self.topo = topo

    def lnLike(self, params): 
        """
        return the log of the likelihood function for the given parameters
        """
        self.setModel(params)
        mod = self.this_model
        lnlike = self.likfun.alikelihood(mod.amplitude, mod.H0,
                                  mod.euler_angles[0],mod.euler_angles[1],
                                  mod.euler_angles[2])
        return lnlike
