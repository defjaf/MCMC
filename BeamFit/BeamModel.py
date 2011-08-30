## maybe should store C^-1 determinant separately??
### get lots of fp errors especially in first iteration of the MCMC


# normalize to 1 at peak
#     could just use inheritance for a version with unit area...

## nb. at and __call__  work with (individual) arrays of x and y!

## functions for getting the beam at a series of params w/ fixed
## location?

## parameters are currently a sequence of sequences, converted from an
## array with 'package' and 'unpackage' functions

## changed: init can now accept 'setup' parameters (maybe should do with metaclass?)
## so may need to setparams() rather than just __init__ in Likelihood

from __future__ import division

from numpy import array, exp, asarray, cos, sin, sqrt, float64
import math
import Proposal


separate_sigma12 = False

if separate_sigma12:
    print "requiring sigma1>sigma2 in prior"


class GaussianBeamModel2D(object):
    """
    model of an unnormalized 2d gaussian beam; parameters
    are location (x, y); and various parameterizations of the shape:
    sigma_major, sigma_minor, angle
    sigma2_x, sigma2_y, rho=corr. coeff.
    """
    
    ### class variables
    centerMin = None
    centerMax = None
    sigMax = None
    
    nparam = 5
    
    fmtstring = "(%.3f %.3f) (%.3f %.3f) %.3f"
    
    #paramBlocks = [0, 0, 1, 1, 2]
    paramBlocks =  [0, 1, 2, 3, 4]
    nBlock = max(paramBlocks)+1
    
    texNames = [r"x", r"y", r"$\sigma_1$", r"$\sigma_2$", r"$\alpha$", r"$A$"]
    
    def setParameters_MajMinAng(self, center, sigmas, angle):
        """set the parameters from x, y, sigma_major, sigma_minor, angle"""
        
        
        self.center = center
        self.sigmas = sigmas
        self.angle = angle % math.pi

        
        c = cos(angle); s = sin(angle)
        s12 = sigmas[0]*sigmas[0]
        s22 = sigmas[1]*sigmas[1]
        
        self.sig2_xy = (c*c*s12 + s*s*s22, s*s*s12 + c*c*s22)
        self.rho = c*s*(s22-s12)/sqrt(self.sig2_xy[0]*self.sig2_xy[1])  # sign issue???
        
        self.set_Cinv()
    
    __init__ = setParameters_MajMinAng  ## default to these params
    
    @classmethod
    def setsigMax(cls, xmax, ymax=None):
        """
        set the maximum sigma (nb. same for x, y)
        """
        
        try:
            cls.sigMax = max(max(xmax.x)-min(xmax.x), max(xmax.y)-min(xmax.y))
        except (NameError, AttributeError):
            if ymax is None:
                ymax = xmax
            cls.sigMax = max(max(xmax)-min(xmax), max(ymax)-min(ymax))
    
    @classmethod
    def setxyRange(cls, xrng, yrng=None, scale=None):
        """
        set the min and max possible values for the x,y coordinates
        nb. this is a classmethod: sets the values for *all* instances!!!
        can set from any of
            a dataset with x and y members
            a sequence's min and max
            a (min, max) tuples
        set the yrange the same as the xrange if the former isn't given
        """
        
        try:
            cls.centerMin = min(xrng.x), min(xrng.y)
            cls.centerMax = max(xrng.x), max(xrng.y)
        except NameError:
            if yrng is None:
                yrng = xrng
            cls.centerMin = (min(xrng), min(yrng))
            cls.centerMax = (max(xrng), max(yrng))
        
        if scale is not None:
            xr = 0.5*(cls.centerMax[0] - cls.centerMin[0])
            yr = 0.5*(cls.centerMax[1] - cls.centerMin[1])
            x0 = 0.5*(cls.centerMax[0] + cls.centerMin[0])
            y0 = 0.5*(cls.centerMax[1] + cls.centerMin[1])
            
            cls.centerMin = x0-scale*xr, y0-scale*yr
            cls.centerMax = x0+scale*xr, y0+scale*yr
            
            
    
    def setParameters_XYRho(self, center, sigma_xy, rho):
        """set the parameters from x, y, sig_x, sig_y, rho=corr.coeff."""
        
        self.center = center
        self.sig2_xy=array(sigma_xy, float64)**2
        self.rho = rho
        self.set_Cinv()
    
    def set_Cinv(self):
        """set the packed inverse correlation matrix values"""
        
        sig2x, sig2y = self.sig2_xy
        det = sig2x * sig2y * (1.0 - self.rho*self.rho)
        
        ## inverse array in a "packed" form (xx, xy, yy)
        self.Cinv = array([
            sig2y, -self.rho*sqrt(sig2x*sig2y), sig2x], float64)/det
    
    def get_XYRho(self):
        """
        get the parameters tuple((x,y), (s2x, s2y), rho)
        nb. for now we return sigma^2, not sigma!
           """
        return self.center, self.sig2_xy, self.rho
    
    def get_MajMinAng(self):
        """ get the parameters tuple((x,y), (s1, s2), angle)
        """
        if self.sigmas is None:
            ## calculate sigma_maj, min; angle from sig_xy, rho
            pass   # for now
        return self.center, self.sigmas, self.angle
    
    def at(self, data):
        """get the [array of] value[s] of the beam for the given dataset"""
        return gauss2d(data.x, data.y, self.center[0], self.center[1], *self.Cinv)
    
    def atxy(self, x, y):
        """
        get the [array of] value[s] of the beam for the given x, y value[s]
        """
        return gauss2d(x, y, self.center[0], self.center[1], *self.Cinv)
    
    __call__ = at
    
    @classmethod
    def prior(cls, center, sigmas, angle):
        """get the unnormalized prior for the parameters
        """
        
        
        if cls.centerMin is not None and (
                center[0] < cls.centerMin[0] or center[1] < cls.centerMin[1] or
                center[0] > cls.centerMax[0] or center[1] > cls.centerMax[1]):
            return 0
        
        if cls.sigMin is None and (sigmas[0]<0 or sigmas[1]<0):
            return 0
        elif sigmas[0]<cls.sigMin or sigmas[1]<cls.sigMin:
            return 0
        
        if cls.sigMax is not None and max(sigmas) > cls.sigMax:
            return 0    # too restrictive?
            
        if separate_sigma12 and sigma2>sigma1:
            return 0
        
        return 1
    
    ## note that we do (angle % pi) in these [probably really only needed in 'package'?]
    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        xy, sig12, ang = param_seqs
        return array( [ xy[0], xy[1], sig12[0], sig12[1], ang % math.pi], dtype=float64)
    
    def package(params_flat):
        """ convert from flat array to structured sequence of parameters """
        return (tuple(params_flat[0:2]), tuple(params_flat[2:4]),
                params_flat[4] % math.pi)
    
    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

## need to do this conversion after we send the methods to the Proposal class
    unpackage=staticmethod(unpackage)
    package=staticmethod(package)

    
    def startfrom(self, data, random=None):
        """
        generate a set of starting parameters for the model:
        non-random version
        center = <x>, <y>
        sigmas = <x^2>-<x>^2, <y^2>-<y>^2
        angle=0
        """
        if random is not None:
            dx = (self.centerMin[0], self.centerMax[0])
            dy = (self.centerMin[1], self.centerMax[1])
            
            start_params = ( (uniform(*dx), uniform(*dy)),
                         (uniform(0,(dx[1]-dx[0])/5), uniform(0,(dy[1]-dy[0])/5)),
                         uniform(0,math.pi/2) )
        else:
            pass
            
            

class GaussianBeamModel2D_xy(GaussianBeamModel2D):
    """like GaussianBeamModel2D, but explicitly use sig_x, sig_y as parameters"""

    texNames = [r"x", r"y", r"$\sigma_x$", r"$\sigma_y$", r"$\rho$", r"$A$"]
    
    def __init__(self, center, sigmas, rho):
        super(GaussianBeamModel2D_xy, self).setParameters_XYRho(center, sigmas, rho)

    @classmethod
    def prior(cls, center, sigmas, rho):
        """get the unnormalized prior for the parameters
        """
        
        if not super(GaussianBeamModel2D_xy, cls).prior(center, sigmas, 0):
            return 0
            
        if rho<-1 or rho>1:
            return 0
            
        return 1
        
    ## note that we do (angle % pi) in these [probably really only needed in 'package'?]
    def unpackage_xy(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        xy, sigxy, rho = param_seqs
        return array( [ xy[0], xy[1], sigxy[0], sigxy[1], rho], dtype=float64)

    def package_xy(params_flat):
        """ convert from flat array to structured sequence of parameters """
        return (tuple(params_flat[0:2]), tuple(params_flat[2:4]), params_flat[4])

    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package_xy,
                                                unpackage=unpackage_xy)

## need to do this conversion after we send the methods to the Proposal class
    unpackage=staticmethod(unpackage_xy)
    package=staticmethod(package_xy)


    def startfrom(self, data, random=None):
        """
        generate a set of starting parameters for the model:
        non-random version
        center = <x>, <y>
        sigmas = <x^2>-<x>^2, <y^2>-<y>^2
        rho=0
        """
        if random is not None:
            dx = (self.centerMin[0], self.centerMax[0])
            dy = (self.centerMin[1], self.centerMax[1])

            start_params = ( (uniform(*dx), uniform(*dy)),
                         (uniform(0,(dx[1]-dx[0])/5), uniform(0,(dy[1]-dy[0])/5)),
                         uniform(-1,1) )
        else:
            pass


def gauss2d(x, y, x0, y0, Cinv_xx, Cinv_xy, Cinv_yy):
    dx = x-x0; dy = y-y0
    return exp(-0.5 * (dx*dx*Cinv_xx + dy*dy*Cinv_yy + 2*dx*dy*Cinv_xy))

        