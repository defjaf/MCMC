## changed: init can now accept 'setup' parameters (maybe should do with metaclass?)
## so may need to setparams() rather than just __init__ in Likelihood

from __future__ import division
from numpy import array, exp, asarray, cos, sin, sqrt, float64
import math
import Proposal
import BeamModel

class GaussianBeamModelRound(BeamModel.GaussianBeamModel2D):
    """
    round (sigma1=sigma2) version of BeamModel.GaussianBeamModel2D
    """

    
    ### class variables
    centerMin = None
    centerMax = None
    sigMax = None

    nparam = 5
    
    fmtstring = "(%.3f %.3f) %.3f"
    
    #paramBlocks = [0, 0, 1]
    paramBlocks =  [0, 1, 2]
    nBlock = max(paramBlocks)+1
    
    texNames = [r"x", r"y", r"$\sigma$"]
    
    def setParameters_MajMin(self, center, sigma):
        """set the parameters from x, y, sigma"""
        
        self.center = center
        self.sigma = sigma

        s2 = sigma**2
        
        self.sig2_xy = (s2, s2)
        self.rho = 0.0
        
        self.set_Cinv()

    __init__ = setParameters_MajMin  ## default to these params
    
    def set_Cinv(self):
        """set the packed inverse correlation matrix values"""

        ## inverse array in a "packed" form (xx, xy, yy)
        self.Cinv = array([1.0/sig2x, 0, 1.0/sig2y], float64)
