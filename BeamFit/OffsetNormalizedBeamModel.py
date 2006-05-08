## should be able to subclass from BeamModel and just add the extra
## parameter as needed and calling the superclass

## nb should define amplitude s.t. \int dx dy modeledBeam(x,y) = Amplitude
##   then can just get stats on Amplitude by marginalizing over all
##   other params


## could really do multiple inheritance with a class for plane (offset, mu, phi)

## New version of 12/2005 which allows offset and 2-d gradient
## NB. MCMC.py & Proposal.py have been modified to allow 'fixed' parameters via
##     sigma=0 proposal density

from __future__ import division
#from numarray import array, exp, asarray, cos, sin, sqrt, float64
from numpy import array, exp, asarray, cos, sin, sqrt, float64, fabs
import math

import Proposal
from BeamModel import GaussianBeamModel2D

class OffsetNormalizedBeamModel(GaussianBeamModel2D):
    """
    model of a *normalized* 2d Gaussian beam model (i.e., it
    integrates to the actual power in data units in the beam)

    now allow for an offset and two-d gradient (i.e. a plane),
    with three parameters:
       offset = value of plane at center
       theta=arccos(mu), phi = sph'l polar coords of normal to the plane
    """
    
    fmtstring = GaussianBeamModel2D.fmtstring + " %.3f (%.3f %.3f %.3f)"
    
    #locBlk = [0, 1, 2, 2]
    locBlk = [0, 1, 2, 3]
    paramBlocks = (GaussianBeamModel2D.paramBlocks 
                   + [GaussianBeamModel2D.nBlock+l for l in locBlk])
    nBlock = max(paramBlocks)+1
    
    texNames = GaussianBeamModel2D.texNames + [r"$O$", r"$\mu$", r"$\phi$"]
    
    nparam = 9

    def __init__(self, center, sigmas, angle, amplitude, plane):
        """
        set the parameters from
            (x, y), (sigma_major, sigma_minor), angle, amplitude,
            (offset, mu, phi)
        """
        self.setParameters_MajMinAng(center, sigmas, angle)
        self.amplitude = amplitude
        sqrtdet = sigmas[0]*sigmas[1]    
        #self.norm = amplitude/sqrtdet/2/math.pi   ### CHECK THIS
        self.norm = amplitude   ### CHECK THIS

        self.offset = plane[0]
        self.offset_mu = plane[1]
        self.offset_phi = plane[2]

    def at(self, data):
        return self.norm*(GaussianBeamModel2D.at(self, data) +
                          self.offset +
                          planexy(self.offset_mu, self.offset_phi, 
                                  data.x-self.center[0], data.y-self.center[1]))

    def atxy(self, x, y):
        return self.norm*(GaussianBeamModel2D.atxy(self, x, y) +
                          self.offset +
                          planexy(self.offset_mu, self.offset_phi, 
                                  x-self.center[0], y-self.center[1]))

    __call__ = at

    @classmethod
    def prior(cls, center, sigmas, angle, amplitude, plane):
        # plane[1]==mu -> only care about 0<theta<pi/2
      #  if amplitude <= 0 or plane[1]<-1 or plane[1]>1:
        if plane[1]<-1 or plane[1]>1:
            return 0
        else:
            ## this calls the superclass method, but accesses the subclass 
            ## class variables
            return super(OffsetNormalizedBeamModel, cls).prior(center, sigmas, angle)


    ## note that we do (angle % pi) in these [probably really only
    ## needed in 'package'?]
        
    ##could probably just append to the
    ## superclass methods' output, but not worth the effort
        
    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        xy, sig12, ang, amp, plane = param_seqs
        return array( [ xy[0], xy[1], sig12[0], sig12[1], ang % math.pi,
                        amp, plane[0], fabs(plane[1]), plane[2] % (2*math.pi)] )
    
    def package(params_flat):
        """ convert from flat array to structured sequence of parameters """
        return (tuple(params_flat[0:2]), tuple(params_flat[2:4]),
                params_flat[4] % math.pi, params_flat[5],
                (params_flat[6], fabs(params_flat[7]), params_flat[8]%(2*math.pi)))

    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

## need to do this conversion after we send the methods to the Proposal class
    unpackage=staticmethod(unpackage)  
    package=staticmethod(package)



### maybe have the 'center'/'origin' as a separate param in theses???
def planexy(mu, phi, x, y):
    """
    return z coord of a plane with normal vector (mu=cos(theta), phi) at (x,y)
    """
    if mu==1: return 0

    return sqrt(1/mu/mu - 1)*(x*cos(phi) + y*sin(phi))



def plane(mu, phi, data):
    """
    return z coord of a plane with normal vector (mu=cos(theta), phi)
    at the data locations
    Not used!
    """

    ## redo the planexy logic here -- faster but ugly
    if mu==1: return 0

    return sqrt(1-mu*mu)/mu*(data.x*cos(phi) + data.y*sin(phi))

