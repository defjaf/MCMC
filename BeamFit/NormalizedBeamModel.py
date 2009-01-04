## should be able to subclass from BeamModel and just add the extra
## parameter as needed and calling the superclass

## nb should define amplitude s.t. \int dx dy modeledBeam(x,y) = Amplitude
##   then can just get stats on Amplitude by marginalizing over all
##   other params


from __future__ import division

from numpy import array, exp, asarray, cos, sin, sqrt, float64
import math

import Proposal

use_xy = False ## True doesn't work yet???

if use_xy:
    from BeamModel import GaussianBeamModel2D_xy as G2D
    print "Using sigma_x, sigma_y, rho params"
else:
    from BeamModel import GaussianBeamModel2D as G2D
    print "Using, sigma_1, sigma_2, angle params"

class NormalizedBeamModel(G2D):
    """
    model of a *normalized* 2d Gaussian beam model (i.e., it
    integrates to the actual power in data units in the beam)
    """
    
    fmtstring = G2D.fmtstring + " %.3f"
    
    nparam = 6

    #print "NB. allowing amplitude < 0!"

    def __init__(self, center, sigmas, angle, amplitude):
        """
        set the parameters from
            (x, y), (sigma_major, sigma_minor), angle, amplitude
        """
        super(OffsetNormalizedBeamModel, self).__init__(center, sigmas, angle)
        self.amplitude = amplitude
        sqrtdet = sigmas[0]*sigmas[1]    
        #self.norm = amplitude/sqrtdet/2/math.pi   ### CHECK THIS
        self.norm = amplitude   ### CHECK THIS

    def at(self, data):
        return self.norm*G2D.at(self, data)

    def atxy(self, x, y):
        return self.norm*G2D.atxy(self, x, y)

    __call__ = at

    @classmethod
    def prior(cls, center, sigmas, angle, amplitude):
#        if amplitude <= 0:
#            return 0
#        else:
#            ## this calls the superclass method, but accesses the subclass 
#            ## class variables
            return super(NormalizedBeamModel, cls).prior(center, sigmas, angle)

    ## note that we do (angle % pi) in these [probably really only needed in 'package'?]
    ## could probably just append to the superclass methods' output, but not worth the effort
    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        xy, sig12, ang, amp = param_seqs
        return array( [ xy[0], xy[1], sig12[0], sig12[1], ang % math.pi, amp] )
    
    def package(params_flat):
        """ convert from flat array to structured sequence of parameters """
        return (tuple(params_flat[0:2]), tuple(params_flat[2:4]),
                params_flat[4] % math.pi, params_flat[5])

    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

## need to do this conversion after we send the methods to the Proposal class
    unpackage=staticmethod(unpackage)  
    package=staticmethod(package)
