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
#from numarray import array, exp, asarray, cos, sin, sqrt, float64
from numpy import array, exp, asarray, cos, sin, sqrt, float64
import math
import Proposal

class ComplexBeamModel(object):
    """
    model of  2d gaussian beam; 
    """
    
    ### class variables

    nparam = 5
    
    fmtstring = "(%.3f %.3f) (%.3f %.3f) %.3f"
    
    #paramBlocks = [0, 0, 1, 1, 2]
    paramBlocks =  range(nparam)
    nBlock = max(paramBlocks)+1
    
    texNames = []
    
    def __init__(self, center, shape):
        self.center = center
        self.params = asarray(shape)

    @classmethod
    def set_numParams(cls, nparam):
        cls.nparam = nparam
        pass
        
    def at(self, data):
        """get the [array of] value[s] of the beam for the given dataset"""
        return beammodel(data.x, data.y, self.params)
        
    def atxy(self, x, y):
        """get the [array of] value[s] of the beam for the given x, y value[s] """
        return beammodel(x, y, self.params)
        
    __call__ = at

    @classmethod
    def prior(cls, center, sigmas, angle):
        """get the unnormalized prior for the parameters
           nb. needs to be staticmethod since classmethod would find the inherited
           instances of the class variables -- which don't propagate
        """
        return 1

    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        return array(param_seqs)
    
    def package(params_flat):
        """ convert from flat array to structured sequence of parameters """
        return array(params_flat)

        
    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

## need to do this conversion after we send the methods to the Proposal class
    unpackage=staticmethod(unpackage)  
    package=staticmethod(package)


def beammodel(x, y, x0, y0, shapeparams):
    pass
