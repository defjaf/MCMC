
from __future__ import division

import math

import matplotlib.pyplot as plt
from numpy import array, exp, asarray, cos, sin, sqrt, float64, linspace, log


import Proposal

## careful -- this is the model for a *single* SED, but [probably]
##            need to marginalize over the ratio r12 for *each object* 
##   should be explicit about r12 for each object 
##    -- make it an array and make this the model for all objects together

### TODO: __call__ for model.quadform (and possibly set_model_vals?)
####   returns matrix of amplitudes for each component at each frequency

h_over_k = 4.799237e-11 ### s K

def blackbody(T, nu):
    """return the blackbody flux at frequency nu for temperature T [CHECK UNITS]"""
    x = h_over_k*nu/T
    prefac = 1.0 #### FIX
    return prefac*nu**3/(exp(x)-1)
    
    
class submmModel2(object):
    """model a submm SED as a two-component grey body: flux = A1 nu^b1 B_nu(T1) + A2 nu^b2 B_nu(T2)
       marginalize analytically over both amplitudes, A1, A2

       needs a __call__ method to return the *matrix* flux(i, nu) for i=1,2!!!!
       and this needs to be correctly handled by the likelihood
       so: __call__ must return in a form usable by likelihood.set_model_vals and ****model.quadform****
         nb. model.quadform takes matrices, A & B, of amplitudes as input and returns A^T N^{-1} B

    """

    nparam = 4
    fmtstring = "%.3f "*4
    paramBlocks =  range(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1

    def __init__(self, b1, T1, b2, T2):

        self.b1 = b1
        self.T1 = T1
        self.b2 = b2
        self.T2 = T2


    @classmethod
    def prior(cls, b1, T1, b2, T2):
        """get the unnormalized prior for the parameters
        """

        if T1<0 or T2<0:
            return 0

        return 1.0

    def at_nu(self, nu):
        return asarray([nu**self.b1 * blackbody(self.T1, nu),
                        nu**self.b2 * blackbody(self.T2, nu)])

    def at(self, data):
        return asarray([data.freq**self.b1 * blackbody(self.T1, data.freq), 
                        data.freq**self.b2 * blackbody(self.T2, data.freq) ]).transpose()
        ### shape=(2,n_freq); is this right? or transpose this?
        
    __call__ = at

    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        return asarray( param_seqs, dtype=float64)


    package = asarray

    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

## need to do this conversion after we send the methods to the Proposal class
    unpackage=staticmethod(unpackage)
    package=staticmethod(package)

    @classmethod
    def startfrom(cls, data=None, random=None):
        """
        generate a set of starting parameters for the model:
        """
        if random is not None:
            cls.start_params = (2., 30., 2., 20.)  ## careful of units
        else:
            pass

        return cls.start_params

    

class submmModel1(object):
    """model a submm SED as a two-component grey body: flux = A1 nu^b1 B_nu(T1) + A2 nu^b2 B_nu(T2)
       marginalize analytically over the overall amplitude, so write as
       unnorm_flux = A (nu^b1 B_nu(T1) + r12 nu^b2 B_nu(T2) )
       logarithmic prior on r12 is symmetric
       
       should just use this as the marginalized one-component model by forcing r12=0?
    """
    
    nparam = 5
    fmtstring = "%.3f "*5
    paramBlocks =  range(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    
    def __init__(self, b1, T1, b2, T2, r12):
        
        self.b1 = b1
        self.T1 = T1
        self.b2 = b2
        self.T2 = T2
        self.r12 = r12
        
        
    def at_nu(self, nu):
        return nu**self.b1 * blackbody(self.T1, nu) + \
               self.r12 * nu**self.b2 * blackbody(self.T2, nu) 
        
    def at(self, data):
        return data.freq**self.b1 * blackbody(self.T1, data.freq) + \
               self.r12 * data.freq**self.b2 * blackbody(self.T2, data.freq) 
        
    __call__ = at    
        
    @classmethod
    def prior(cls, b1, T1, b2, T2, r12):
        """get the unnormalized prior for the parameters
        """

        if T1<0 or T2<0:
            return 0

        return log(r12)
        
        

    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        return asarray( param_seqs, dtype=float64)

    package = tuple

    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

## need to do this conversion after we send the methods to the Proposal class
    unpackage=staticmethod(unpackage)
    package=staticmethod(package)

    def plot(self, data):
        """plot the data and the model"""
        f = np.linspace(data.freq[0], data.freq[-1], 100)
        model_flux = self.at(f)
        data.plot()
        plt.plot(f, model_flux)
        
        
    @classmethod
    def startfrom(cls, data, random=None):
        """
        generate a set of starting parameters for the model:
        """
        if random is not None:
            start_params = (2., 10., 2., 5., 1.0)  ## careful of units
        else:
            pass
            
            
class submmModel_normalized(submmModel1):
    """model a submm SED as a two-component grey body: flux = A1 nu^b1 B_nu(T1) + A2 nu^b2 B_nu(T2)

    """

    nparam = 6
    fmtstring = "%.3f "*6
    paramBlocks =  range(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1

    def __init__(self, A1, b1, T1, A2, b2, T2):

        self.A1 = A1
        self.b1 = b1
        self.T1 = T1
        self.b2 = b2
        self.T2 = T2
        self.A2 = A2


    def at_nu(self, nu):
        return self.A1 * nu**self.b1 * blackbody(self.T1, nu) + \
               self.A2 * nu**self.b2 * blackbody(self.T2, nu) 

    def at(self, data):
        return self.A1 * data.freq**self.b1 * blackbody(self.T1, data.freq) + \
               self.A2 * data.freq**self.b2 * blackbody(self.T2, data.freq) 

    __call__ = at    

    @classmethod
    def prior(cls, A1, b1, T1, A2, b2, T2):
        """get the unnormalized prior for the parameters
        """

        if T1<0 or T2<0:
            return 0

        return 1

    def plot(self, data):
        """plot the data and the model"""
        f = np.linspace(data.freq[0], data.freq[-1], 100)
        model_flux = self.at(f)
        data.plot()
        plt.plot(f, model_flux)



    def startfrom(self, data, random=None):
        """
        generate a set of starting parameters for the model:
        """
        if random is not None:
            start_params = (1., 2., 10., 1., 2., 5.)  ## careful of units
        else:
            pass

