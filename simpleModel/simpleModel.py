""" """

## make this the base class????

## nb. really 2 kinds of models:
##  1. the parametrized model depends on the data (i.e., you need to
##     find at a specific seto of points depening on where you've
##     observed)
##  2. the parameterized model is 'global', calculated once
##     and for all (e.g., the CMB case, although this is only because the
##     likelihood function works directly on the C_l and carries around
##     the dataset information)


from __future__ import division

import math
from operator import isNumberType, __and__,  isSequenceType
from copy import copy

from numarray import array, arange, Float64, zeros
import Proposal

class simpleModel(object):
    """ general polynomial model for the data;
        params are the list [a_i,...] in y = a0 + a1 x + a2 x^2 + ... + noise
        nb. need to switch to scipy.poly1d object...
    """

    def __init__(self, params):
        self.params=params
    
    def __call__(self, data):
        return self.atx(data.x)

    def atx(self, xarr):
        """ return the parametrized model for these params
            currently a really terrible numerical implementation
        """
        ## can probably do this in a single sum using an axis= parameter?
        return array([ sum( a*x**n for (n, a) in enumerate(self.params) )
                       for x in xarr ])

    @staticmethod
    def prior(*params):
        return 1.0

    def package(par): return par
    def unpackage(par): return par
    
    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

    unpackage=staticmethod(unpackage)  
    package=staticmethod(package)
    
