



import math

import matplotlib.pyplot as plt
from numpy import (array, exp, asarray, cos, sin, sqrt, float64, linspace, log, errstate, 
                    min, max, vectorize, expm1, zeros)
from scipy import special,integrate
#import numexpr as ne

speed_of_light = 299792.458 ### micron GHz

singleParams = True  ## to vary each parameter individually
if not singleParams: print("varying parameters together!")

from .. import Proposal

## careful -- this is the model for a *single* SED, but
##            would need to marginalize over the ratio r12 for *each object* 
##   should be explicit about r12 for each object 
##    -- make it an array and make this the model for all objects together


#### for unnormed cases, amplitudes are calculated in the associated likelihood class
####  (not checked for more than one amplitude)
####  TODO: use amplitude>0 prior???

##### TODO: get the total far-IR luminosity [for each component] of an object with a given model DONE
##### NOTE: the integral is done dnu in GHz so need to multiply by 1e9 to convert to (Jy Hz)

##### DONE: only plot data once (the first time?) on a single plot (but it's done by hand)
##### TODO: check that this is correct in driver.py?

##### TODO: correctly normalized prior (done?)



#### AHJ 08/2015 -- for two-component models, try to reorder the parameters rather than enforce via prior
#####                do this at the "__init__" step -- not sure this is right...*****
#####               need to do this in "package"? 
#####                  (See e.g., OffsetNormalizedBeamModel which does a "% pi" operation there)
#####                  DONE

#####               rewrite at(self, data) to use the data to normalize the amplitudes 
#####                   model = A * data(nu_bar) * F(nu)/F(nu_bar)
#####                   should be possible since it has access to the data
#### AHJ 09/2015 -- do normalization as above for optical depth model
#####                   but now, when nu_0 -> infty, P(nu_0) -> const.
#####               rescale nu_0 -> nu_0/1000 so O(1) 
#####            both of these done for pystan runs -- not checked in MCMC runs
                    

h_over_k = 0.04799237 ###  K/Ghz
prefac = 1.0e-9 #### FIXME: find a physical definition to go here
nu_b = 1000.0   ### frequency at which to normalize the SEDs

minTemp, maxTemp = 3.0, 100.0
print("min Temp = %f K; max Temp = %f K" % (minTemp, maxTemp))
minb, maxb = 0., 3.

def startfrom_generic(start, stds, posidx=(), random=True):
    return start
# randomization now handled elsewhere.
#     while True:
#         if random:
#             start += np.random(len(stds))*stds
#         if all(start[p]>0 for p in posidx):
#             return start

# try:
# 
#     from blackbody import greybody, blackbody, lux
# 
#     print("Got blackbody.pyx")
# 
# except ImportError:
# print("python version of blackbody")
def blackbody(T, nu, nu_norm=False):
    """return the blackbody flux at frequency nu for temperature T [CHECK UNITS]"""

    # if T==0:
    #     return 0
    
    x = h_over_k*nu/T
    
    with errstate(over='ignore'):
      #  return ne.evaluate("prefac*nu**3/(exp(x)-1)")
        if nu_norm:
            x_norm = h_over_k*nu_norm/T
            return (nu/nu_norm)**3 * expm1(x_norm)/expm1(x)
        else:
            return prefac * nu**3/expm1(x)
      
      
def greybody(beta, T, nu, nu_norm=False):
    """return the greykody flux at frequency nu for temperature T [CHECK UNITS]"""

    # if T==0:
    #     return 0

    x = h_over_k*nu/T
    with errstate(over='ignore'):
      #  return ne.evaluate("prefac*nu**3/(exp(x)-1)")
#          return prefac*nu**(3+beta)/(exp(x)-1)
        if nu_norm:
            x_norm = h_over_k*nu_norm/T
            return (nu/nu_norm)**(3+beta) * expm1(x_norm)/expm1(x)
        else:
            return prefac * nu_b**(-beta) * nu**(3+beta)/expm1(x)

@vectorize
def totalflux(beta, T, nu1=None, nu2=None, nu_norm=False):
  """
  calculate the total flux of a grey body (with prefactors defined as above) over (nu1,nu2)
  NOTE: the integral is done dnu in GHz so need to multiply by 1e9 to convert to (Jy Hz)
  THIS IS INCORRECT FOR self.amplitude DEFINED WITH nu_norm
  """

  if nu1 is None and nu2 is None:
      ## return analytic expression for bolometric flux
      return prefac * nu_b**(-beta) * (T/h_over_k)**(4+beta) * \
             special.gamma(4+beta)*special.zeta(4+beta,1)
  else:
      ## do numeric integral
      return integrate.quad(lambda nu: greybody(beta, T, nu, nu_norm=nu_norm), nu1, nu2)[0]

    
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
    paramBlocks =  list(range(nparam)) if singleParams else zeros(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$\beta_1$", r"$T_1$", r"$\beta_2$", r"$T_2$"]

    def __init__(self, b1, T1, b2, T2):
    
        ### AHJ 08/2015 keep T1<T2 always
#         if T1>T2:
#             T1, T2 = T2, T1
#             b1, b2 = b2, b1    

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
            
        ## want to separate the two cases: force T1<T2
#         if T1>T2: return 0   # AHJ 08/2015

        return 1.0

    def at_nu(self, nu):
        return asarray([greybody(self.b1, self.T1, nu),
                        greybody(self.b2, self.T2, nu)])

    def at(self, data):
        return asarray([greybody(self.b1, self.T1, data.freq), 
                        greybody(self.b2, self.T2, data.freq) ]).transpose()
        
    __call__ = at

    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        return asarray( param_seqs, dtype=float64)


#     package = asarray
    def package(params_flat):
        b1, T1, b2, T2 = params_flat
        if T1>T2:
            params_flat = array( (b2, T2, b1, T1) )
        return asarray(params_flat)

    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

## need to do this conversion after we send the methods to the Proposal class
    unpackage=staticmethod(unpackage)
    package=staticmethod(package)

    @classmethod
    def startfrom(cls, data=None, random=None):
        """
        generate a set of starting parameters for the model: b1, T1, b2, T2
        """
        start_params = (2., 20., 2., 30.)  ## careful of units
        cls.start_stds = stds = (0.5, 6.0, 0.5, 4.0)
        posidx = (1,3)
        cls.start_params = startfrom_generic(start_params, stds, posidx, random=random)
        return cls.start_params


class submmModel1(object):
    """model a submm SED as a one-component grey body: flux = A nu^b B_nu(T)
       marginalize analytically over the overall amplitude
    """

    nparam = 2
    fmtstring = "%.3f "*2
    paramBlocks =  list(range(nparam)) if singleParams else zeros(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$\beta$", "$T$"]

    def __init__(self, b, T):

        self.b = b
        self.T = T


    def at_nu(self, nu):
        return greybody(self.b, self.T, nu) 

    def at(self, data):
        return greybody(self.b, self.T, data.freq) 
        
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
        start_params = (2.0, 10.0)
        cls.start_stds = stds = (0.5, 3.0)
        posidx = (1,)
        cls.start_params = startfrom_generic(start_params, stds, posidx, random=random)
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
    paramBlocks =  list(range(nparam)) if singleParams else zeros(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    
    def __init__(self, b1, T1, b2, T2, r12):
        
        self.b1 = b1
        self.T1 = T1
        self.b2 = b2
        self.T2 = T2
        self.r12 = r12
        
        
    def at_nu(self, nu):
        return greybody(self.b1, self.T1, nu) + \
               self.r12 * greybody(self.b2, self.T2, nu) 
        
    def at(self, data):
        return greybody(self.b1, self.T1, data.freq) + \
               self.r12 * greybody(self.b2, self.T2, data.freq) 
        
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

    def plot(self, data, logplot=True, label=None, plot_data=True):
        """plot the data and the model"""
        f = linspace(min(data.freq), max(data.freq), 100)
        model_flux = self.at_nu(f)
        if plot_data:
            data.plot(logplot=logplot)
        plt.plot(f, model_flux, label=label)
        
        
    @classmethod
    def startfrom(cls, data=None, random=None):
        """
        generate a set of starting parameters for the model: b1, T1, b2, T2, r12
        """
        start_params = (2., 5., 2., 10., 1.0)  ## careful of units
        cls.start_stds = stds = (0.5, 3, 0.5, 2.0, 0.25)
        posidx = (1,3,4)
        cls.start_params = startfrom_generic(start_params, stds, posidx, random=random)
        return cls.start_params

            
class submmModel2_normalized(object):
    """model a submm SED as a two-component grey body: flux = A1 nu^b1 B_nu(T1) + A2 nu^b2 B_nu(T2)

    """

    nparam = 6
    fmtstring = "%.3f "*6
    paramBlocks =  list(range(nparam)) if singleParams else zeros(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$A_1$", r"$\beta_1$", r"$T_1$", r"$A_2$", r"$\beta_2$", r"$T_2$"]

    def __init__(self, A1, b1, T1, A2, b2, T2):

        ### AHJ 08/2015 keep T1<T2 always
#         if T1>T2:
#             T1, T2 = T2, T1
#             b1, b2 = b2, b1
#             A1, A2 = A2, A1

        self.A1 = A1
        self.b1 = b1
        self.T1 = T1
        self.b2 = b2
        self.T2 = T2
        self.A2 = A2
        

    def at_nu(self, nu):
        return self.A1 * greybody(self.b1, self.T1, nu, nu_norm=nu_b) + \
               self.A2 * greybody(self.b2, self.T2, nu, nu_norm=nu_b) 

    def at(self, data):
        return self.A1 * greybody(self.b1, self.T1, data.freq, nu_norm=nu_b) + \
               self.A2 * greybody(self.b2, self.T2, data.freq, nu_norm=nu_b) 

    __call__ = at    

    
    def flux(self, nu1=None, nu2=None):
        """return the flux between the given frequencies for each component"""
        return asarray([self.A1*totalflux(self.b1, self.T1, nu1, nu2, nu_norm=nu_b), 
                        self.A2*totalflux(self.b2, self.T2, nu1, nu2, nu_norm=nu_b)])
        

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
#         if T1>T2: return 0   # AHJ 08/2015 removed
        

        return 1

    def plot(self, data, wavelength=True, logplot=True, label=None, plot_data=True):
        """plot the data and the model"""
        f = linspace(min(data.freq), max(data.freq), 100)
        model_flux = self.at_nu(f)
        if plot_data:
            data.plot(fmt='o', wavelength=wavelength, logplot=logplot)
        if wavelength:
            f = speed_of_light/f
        plt.plot(f, model_flux, label=label)

#     package = asarray
# 
    def package(params_flat):
        A1, b1, T1, A2, b2, T2 = params_flat
        if T1>T2:
            params_flat = asarray( (A2, b2, T2, A1, b1, T1) )
        return asarray(params_flat)

        

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
        generate a set of starting parameters for the model: A1, b1, T1, A2, b2, T2
        """
        start_params = (1., 2., 5., 1., 2., 10.)  ## careful of units
        cls.start_stds = stds = (0.25, 0.5, 5.0, 0.25, 0.5, 5.0)
        posidx = (0,2,3,5)
        cls.start_params = startfrom_generic(start_params, stds, posidx, random=random)
        return cls.start_params            

class submmModel2_normalized_logA(submmModel2_normalized):
    """model a submm SED as a two-component grey body: flux = 10**logA1 nu^b1 B_nu(T1) + 10**A2 nu^b2 B_nu(T2)
    """
    texNames = [r"$\log A_1$", r"$\beta_1$", r"$T_1$", r"$\log A_2$", r"$\beta_2$", r"$T_2$"]

    def __init__(self, logA1, b1, T1, logA2, b2, T2):
    
        ### AHJ 08/2015 keep T1<T2 always
#         if T1>T2:
#             T1, T2 = T2, T1
#             b1, b2 = b2, b1
#             logA1, logA2 = logA2, logA1
# 

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
#         if T1>T2: return 0  # AHJ 08/2015
        
        return 1

    def package(params_flat):
        logA1, b1, T1, logA2, b2, T2 = params_flat
        if T1>T2:
            params_flat = asarray( (logA2, b2, T2, logA1, b1, T1) )
        return asarray(params_flat)

    def unpackage(param_seqs):
        """ convert from structured sequence of parameters to flat array """
        return asarray( param_seqs, dtype=float64)

    ## nb. an *instance* of proposal; should pass the class [name] to this?
    proposal = Proposal.GenericGaussianProposal(package=package,
                                                unpackage=unpackage)

## need to do this conversion after we send the methods to the Proposal class
    unpackage=staticmethod(unpackage)
    package=staticmethod(package)


class submmModel1_normalized(submmModel2_normalized):
    """model a submm SED as a one-component grey body: flux = A nu^b B_nu(T)

    """

    nparam = 3
    fmtstring = "%.3f "*3
    paramBlocks =  list(range(nparam)) if singleParams else zeros(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$A$", r"$\beta$", r"$T$"]
    
    def __init__(self, A, b, T):

        self.A = A
        self.b = b
        self.T = T


    def at_nu(self, nu):
        return self.A * greybody(self.b, self.T, nu, nu_norm=nu_b)

    def at(self, data):
        return self.A * greybody(self.b, self.T, data.freq, nu_norm=nu_b)

    __call__ = at    

    def flux(self, nu1=None, nu2=None):
        """return the flux between the given frequencies for each component"""
        return array([self.A*totalflux(self.b, self.T, nu1, nu2)])   ### return a scalar?


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
        start_params = (1., 2., 10.)  ## careful of units
        cls.start_stds = stds = (0.25, 0.5, 3.0)
        posidx = (0,2)
        cls.start_params = startfrom_generic(start_params, stds, posidx, random=random)
        return cls.start_params
        
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


        
class submmModel1_opticallythick(submmModel1_normalized):
    """model a submm SED as an optically-thick black body: flux = A (1-exp(-tau)) B_nu(T)
            tau = (nu/nu_0)^b
            
            normalize amplitude to nu=nu_bar
            rescale nu_0 to nu_0/1000.0 (i.e., THz?)
                done for pystan, and still works for the MCMC code
               
            makes P(nu_0) ~ const as nu_0 -> infty
                need prior on nu_0 ~ exponential(1/3)?  or Normal(1,3)?
                
    """

    nparam = 4   # A, b, T, nu_0
    fmtstring = "%.3f "*nparam
    paramBlocks = list(range(nparam)) if singleParams else zeros(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$A$", r"$\beta$", r"$T$", r"$\nu_0$"]
    nu_rescale = 1000.0   ### rescale so O(1)
    nu_renorm = True      ### normalize to nu=nu_b 
    
    def __init__(self, A, b, T, nu_0):

        self.A = A
        self.b = b
        self.T = T
        self.nu_0 = nu_0

    if nu_renorm:
        def at_nu(self, nu):
            tau = (nu/self.nu_0/self.nu_rescale)**self.b
            tau_b = (nu_b/self.nu_0/self.nu_rescale)**self.b
            return self.A * (1.0-exp(-tau))/(1.0-exp(-tau_b)) * blackbody(self.T, nu, nu_norm=nu_b)

        def at(self, data):
            tau = (data.freq/self.nu_0/self.nu_rescale)**self.b
            tau_b = (nu_b/self.nu_0/self.nu_rescale)**self.b
            return self.A * (1.0-exp(-tau))/(1.0-exp(-tau_b)) * blackbody(self.T, data.freq, nu_norm=nu_b)
    else:
        def at_nu(self, nu):
            tau = (nu/self.nu_0/nu_rescale)**self.b
            return self.A * (1.0-exp(-tau)) * blackbody(self.T, nu, nu_norm=nu_b)

        def at(self, data):
            tau = (data.freq/self.nu_0/1000.0)**self.b
            return self.A * (1.0-exp(-tau)) * blackbody(self.T, data.freq, nu_norm=nu_b)

    __call__ = at    

    def flux(self, nu1=None, nu2=None):
        """return the flux between the given frequencies for each component"""
        return array([self.A*totalflux(self.b, self.T, nu1, nu2, nu_norm=nu_b)])   ### return a scalar?


    @classmethod
    def prior(cls, A, b, T, nu_0):
        """get the unnormalized prior for the parameters"""
        
        if A<0: 
            return 0

        if T<minTemp or T>maxTemp:
            return 0
            
        if b<minb or b>maxb:
            return 0
            
        if nu_0<0:
            return 0
        elif cls.nu_rescale and cls.nu_renorm:
            #return exp(-0.5*(nu_0-1.0)/3.0**2)
            return exp(-3.0*nu_0)

        return 1

        
    @classmethod
    def startfrom(cls, data=None, random=None):
        """
        generate a set of starting parameters for the model:
        """
        start_params = (1., 2., 10., 1.0)  ## careful of units
        cls.start_stds = stds = (0.25, 0.5, 3.0, 0.3)
        posidx = (0,2,3)
        cls.start_params = startfrom_generic(start_params, stds, posidx, random=random)
        return cls.start_params


        
class submmModel1_opticallythick_logA(submmModel1_normalized):
    """model a submm SED as an optically-thick black body: flux = A (1-exp(-tau)) B_nu(T)
            tau = (nu/nu_b)^b

    """

    nparam = 4   # A, b, T, [nu_0?]
    fmtstring = "%.3f "*nparam
    paramBlocks =  list(range(nparam)) if singleParams else zeros(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"log $A$", r"$\beta$", r"$T$", r"$\nu_0$"]
    nu_rescale = 1000.0   ### rescale so O(1)
    nu_renorm = True      ### normalize to nu=nu_b 
    
    def __init__(self, logA, b, T, nu_0):

        self.A = 10.0**logA
        self.b = b
        self.T = T
        self.nu_0 = nu_0

    @classmethod
    def prior(cls, logA, b, T, nu_0):
        """get the unnormalized prior for the parameters"""
        
        if T<minTemp or T>maxTemp:
            return 0
            
        if b<minb or b>maxb:
            return 0
            
        if nu_0<0:
            return 0
        elif cls.nu_rescale and cls.nu_renorm:
            #return exp(-0.5*(nu_0-1.0)/3.0**2)
            return exp(-3.0*nu_0)

        return 1        

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

