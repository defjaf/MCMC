
from __future__ import division

import math

import matplotlib.pyplot as plt
from numpy import array, exp, asarray, cos, sin, sqrt, float64, linspace, log, errstate, min, max

#import numexpr as ne

speed_of_light = 299792.458 ### micron GHz

import Proposal

## careful -- this is the model for a *single* SED, but
##            would need to marginalize over the ratio r12 for *each object* 
##   should be explicit about r12 for each object 
##    -- make it an array and make this the model for all objects together


#### for unnormed cases, amplitudes are calculated in the associated likelihood class
####  (not checked for more than one amplitude)
####  TODO: use amplitude>0 prior???

h_over_k = 0.04799237 ###  K/Ghz
prefac = 1.0e-10 #### FIXME: find a physical definition to go here

minTemp, maxTemp = 3.0, 100.0
print "min Temp = %f K; max Temp = %f K" % (minTemp, maxTemp)
minb, maxb = 0., 3.

def blackbody(T, nu):
    """return the blackbody flux at frequency nu for temperature T [CHECK UNITS]"""
    
    # if T==0:
    #     return 0
        
    x = h_over_k*nu/T
    with errstate(over='ignore'):
      #  return ne.evaluate("prefac*nu**3/(exp(x)-1)")
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
    texNames = [r"$\beta_1$", r"$T_1$", r"$\beta_2$", r"$T_2$"]

    def __init__(self, b1, T1, b2, T2):

        self.b1 = b1
        self.T1 = T1
        self.b2 = b2
        self.T2 = T2


    @classmethod
    def prior(cls, b1, T1, b2, T2):
        """get the unnormalized prior for the parameters
        """

        if T1<minTemp or T2<minTemp:
            return 0
            
        if b1<minb or b2<minb or b1>maxb or b2>maxb:
            return 0
            
        ### want to separate the two cases: force T1<T2
        if T1>T2: return 0

        return 1.0

    def at_nu(self, nu):
        return asarray([nu**self.b1 * blackbody(self.T1, nu),
                        nu**self.b2 * blackbody(self.T2, nu)])

    def at(self, data):
        return asarray([data.freq**self.b1 * blackbody(self.T1, data.freq), 
                        data.freq**self.b2 * blackbody(self.T2, data.freq) ]).transpose()
        
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
    """model a submm SED as a one-component grey body: flux = A nu^b B_nu(T)
       marginalize analytically over the overall amplitude
    """

    nparam = 2
    fmtstring = "%.3f "*2
    paramBlocks =  range(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$\beta$", "$T$"]

    def __init__(self, b, T):

        self.b = b
        self.T = T


    def at_nu(self, nu):
        return nu**self.b * blackbody(self.T, nu) 

    def at(self, data):
        return data.freq**self.b * blackbody(self.T, data.freq) 
        
    __call__ = at    

    @classmethod
    def prior(cls, b, T):
        """get the unnormalized prior for the parameters
        """

        if T<minTemp:
            return 0

        if b<minb or b>maxb:
            return 0

        return 1.

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
            cls.start_params = (2.0, 10.0)  ## careful of units
        else:
            pass
        return cls.start_params


class submmModel_ratio(object):
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

        if b1<minb or b2<minb or b1>maxb or b2>maxb:
            return 0

        if T1<minTemp or T2<minTemp:
            return 0
            
        return log(r12)
        

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

    def plot(self, data, logplot=True):
        """plot the data and the model"""
        f = linspace(data.freq[0], data.freq[-1], 100)
        model_flux = self.at_nu(f)
        data.plot(logplot=logplot)
        plt.plot(f, model_flux)
        
        
    @classmethod
    def startfrom(cls, data=None, random=None):
        """
        generate a set of starting parameters for the model:
        """
        if random is not None:
            cls.start_params = (2., 10., 2., 5., 1.0)  ## careful of units
        else:
            pass
        return cls.start_params
            
            
class submmModel2_normalized(object):
    """model a submm SED as a two-component grey body: flux = A1 nu^b1 B_nu(T1) + A2 nu^b2 B_nu(T2)

    """

    nparam = 6
    fmtstring = "%.3f "*6
    paramBlocks =  range(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$A_1$", r"$\beta_1$", r"$T_1$", r"$A_2$", r"$\beta_2$", r"$T_2$"]

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

        if T1<minTemp or T2<minTemp or A1<0 or A2<0:
            return 0
            
        if T1>maxTemp or T2>maxTemp:
            return 0
            
        if b1<minb or b2<minb or b1>maxb or b2>maxb:
            return 0

        ### want to separate the two cases: force T1<T2
        if T1>T2: return 0
        

        return 1

    def plot(self, data, wavelength=True, logplot=True):
        """plot the data and the model"""
        f = linspace(min(data.freq), max(data.freq), 100)
        model_flux = self.at_nu(f)
        data.plot(fmt='o', wavelength=wavelength, logplot=logplot)
        if wavelength:
            f = speed_of_light/f
        plt.plot(f, model_flux)

    package = asarray

    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        return asarray( param_seqs, dtype=float64)

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
            cls.start_params = (1., 2., 10., 1., 2., 5.)  ## careful of units
        else:
            pass
            

class submmModel2_normalized_logA(submmModel2_normalized):
    """model a submm SED as a two-component grey body: flux = 10**logA1 nu^b1 B_nu(T1) + 10**A2 nu^b2 B_nu(T2)
    """
    texNames = [r"$\log A_1$", r"$\beta_1$", r"$T_1$", r"$\log A_2$", r"$\beta_2$", r"$T_2$"]

    def __init__(self, logA1, b1, T1, logA2, b2, T2):

        self.A1 = 10.0**logA1
        self.b1 = b1
        self.T1 = T1
        self.b2 = b2
        self.T2 = T2
        self.A2 = 10.0**logA2

    @classmethod
    def prior(cls, logA1, b1, T1, logA2, b2, T2):
        """get the unnormalized prior for the parameters
        """

        if T1<minTemp or T2<minTemp:
            return 0
            
        if T1>maxTemp or T2>maxTemp:
            return 0
            
        if b1<minb or b2<minb or b1>maxb or b2>maxb:
            return 0

        ### want to separate the two cases: force T1<T2
        if T1>T2: return 0
        
        return 1


class submmModel1_normalized(submmModel2_normalized):
    """model a submm SED as a one-component grey body: flux = A nu^b B_nu(T)

    """

    nparam = 3
    fmtstring = "%.3f "*3
    paramBlocks =  range(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$A$", r"$\beta$", r"$T$"]
    
    def __init__(self, A, b, T):

        self.A = A
        self.b = b
        self.T = T


    def at_nu(self, nu):
        return self.A * nu**self.b * blackbody(self.T, nu)

    def at(self, data):
        return self.A * data.freq**self.b * blackbody(self.T, data.freq)

    __call__ = at    

    @classmethod
    def prior(cls, A, b, T):
        """get the unnormalized prior for the parameters"""

        if T<minTemp or A<0 or T>maxTemp:
            return 0
            
        if b<minb or b>maxb:
            return 0

        return 1
        
        #     def plot(self, data):
        #         """plot the data and the model"""
        #         f = linspace(data.freq[0], data.freq[-1], 100)
        #         model_flux = self.at(f)
        #         data.plot()
        #         plt.plot(f, model_flux)
        # 
        #     package = asarray
        # 
        #     def unpackage(param_seqs):
        #         """ convert from structured sequence of parameters to flat array """
        #         return asarray( param_seqs, dtype=float64)
        # 
        #     ## nb. an *instance* of proposal; should pass the class [name] to this?
        #     proposal = Proposal.GenericGaussianProposal(package=package,
        #                                                 unpackage=unpackage)
        # 
        # ## need to do this conversion after we send the methods to the Proposal class
        #     unpackage=staticmethod(unpackage)
        #     package=staticmethod(package)
        
    @classmethod
    def startfrom(cls, data=None, random=None):
        """
        generate a set of starting parameters for the model:
        """
        if random is None:
            cls.start_params = (1., 2., 10.)  ## careful of units
        else:
            pass

class submmModel1_normalized_logA(submmModel1_normalized):
    """model a submm SED as a one-component grey body: flux = 10**(logA) nu^b B_nu(T)
    """
    texNames = [r"$\log A$", r"$\beta$", r"$T$"]

    def __init__(self, logA, b, T):

        self.A = 10.0**logA
        self.b = b
        self.T = T

    @classmethod
    def prior(cls, logA, b, T):
        """get the unnormalized prior for the parameters"""

        if T<minTemp or T>maxTemp:
            return 0

        if b<minb or b>maxb:
            return 0

        return 1

