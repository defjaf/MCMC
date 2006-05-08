""" represents a topologically non-trivial, slightly-closed universe """

from __future__ import division

import math

from numpy import array, any
import Proposal

class topo_model(object):
    """
    model for a topologically non-trivial (non simply connected) universe
    
    parameters:
        fluctuation amplitude
        euler angles (alpha, beta, gamma)
        hubble constant (only appropriate for these slightly-closed universes
           where H_0 determines the curvature scale and thus the allowed topologies

    """
    
    fmtstring = "%.3f (%.3f %.3f %.3f) %.3f"

    nparam = 5
    
    #paramBlocks = [0, 1, 1, 1, 2]
    paramBlocks = range(nparam)
    nBlock = max(paramBlocks)+1
    
    texNames = [r"$A$", r"$\alpha$", r"$\beta$", r"$\gamma$", r"$H_0$"]
    
    H0_min = 0.0
    H0_max = 100.0
    
    ang_lims = array([2.0, 1.0, 2.0])*math.pi

    def __init__(self, amplitude, euler_angles, H0):
        """
        initialize the parameters
        """
        self.amplitude = amplitude
        self.euler_angles = euler_angles
        self.H0 = H0

    def __call__(self):
        """
        In future may want to return a full correlation matrix here?
        for now, don't do anything!
        """
        return None
            
    @classmethod
    def prior(cls, amplitude, euler_angles, H0):

        if amplitude < 0:
            return 0
        if cls.H0_min is not None and H0<cls.H0_min:
            return 0
        if cls.H0_max is not None and H0>cls.H0_max:
            return 0
            
        if any(euler_angles<0):
            return 0
        if any(euler_angles>cls.ang_lims):
            return 0
        
        return 1

    ## note that we do (angle % pi) in these [probably really only needed in 'package'?]
    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        amp, euler, H0 = param_seqs
        
#        euler %= self.ang_lims   ### wrap around to inside the limits
        
        return array( [ amp, 
                        euler[0], euler[1], euler[2],
                        H0] )
    
    def package(params_flat):
        """ convert from flat array to structured sequence of parameters """
#        return (params_flat[0], tuple(array(params_flat[1:4])%self.ang_lims),
        return (params_flat[0], tuple(params_flat[1:4]),
                params_flat[4] )
                    
    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

    unpackage=staticmethod(unpackage)  
    package=staticmethod(package)
