
from __future__ import division

import math

import matplotlib.pyplot as plt
from numpy import array, exp, asarray, cos, sin, sqrt, float64, linspace, log, errstate, min, max, vectorize, expm1, logspace
from scipy import special,integrate
import numpy as np
#import numexpr as ne

speed_of_light = 299792458 ### meter / Second
c2 = speed_of_light**2
kB = 1.38065e-23 ## Joule/Kelvin  ### boltzmann
h = 6.62607e-34  ## Joule Second ### Planck
Tcmb = 2.72548   # CMB temperature (K)
solid_angle = 0.00382794 # Sr

### how can I parameterize the model rather than hardcode the prior as below?

## nb. 1 Jy = 1e-26 Joule / m^2 
 
import Proposal

from model import greybody, submmModel2_normalized

# 
# The model that I'm currently fitting is (in IDL code):
# 
# Synchrotron (simple power law):
# 
# S = sync_amp * (X)^(sync_spec)
# 
# Free-free (complicated, Te is fixed at 8000K)
# 
# g_ff = alog(4.955d-2 / nu) + 1.5*alog(Te)
# tau_ff = 3.014d-2 * Te^(-1.5) * nu^(-2d) * EM * g_ff
# T_ff = Te * (1d - exp(-tau_ff))
# 
# ; deal with very small tau values in the exponential
# smalltau = where(tau_ff LT 1.0e-10)
# IF (smalltau[0] GE 0) THEN T_ff[smalltau] = Te * (tau_ff[smalltau] - (-tau_ff[smalltau])^2d/ 2d - (-tau_ff[smalltau]^3d)/6d)
# 
# S = 2d * 1381d * T_ff * (nu*1d9)^2 / c^2  * solid_angle
# 
# and the spinning dust model is attached (it needs interpolation).
# 
# 
# There's also CMB, which we're modeling as:
# tcmb = 2.726   ; CMB temperature (K)
# 
# xx = h * X*1.0d9 / (k * tcmb)
# xx2 = h * X*1.0d9 / (k * (tcmb+deltaT))
# S = (2.0e26 * h * (X*1d9)^3 * solid_angle  / c^2) * ((1.0/(exp(xx2)-1d)) - (1.0/(exp(xx)-1d)))
# 
# The solid_angle for this aperture is 0.00382794.

def setup_AME(fname="M31/spdust2_wim.dat"):
    nu_GHz, flux = np.loadtxt(fname, unpack=True)
    return np.log(nu_GHz), np.log(flux)

lognu_AME, logflux_AME = setup_AME()

def AME(nu_GHz, lognu_AME=lognu_AME, logflox_AME=logflux_AME):
    return np.exp(np.interp(np.log(nu_GHz), lognu_AME, logflux_AME))

def freefree(EM, nu_GHz, Te=8000.0, Omega=solid_angle):
    ### nb, [nu] = GHz; output in Jy = 1e26 Joule / m^2
    nu2 = nu_GHz * nu_GHz
    g_ff = np.log(4.955e-2/nu_GHz) + 1.5*np.log(Te)
    tau_ff = 3.014e-2 * Te**(-1.5) * g_ff * EM / nu2
    Tff = -Te*np.expm1(-tau_ff)  # nb. expm1(x) = exp(x)-1 
    Sff = 2*kB*Tff * Omega * nu2*1e18/c2  ## 1e18 converts nu2 to Hz^2
    return 1e26*Sff  ### 1e26 converts to Jy
    
def Bnu(nu_Hz, T_K=Tcmb):
    """ black body, in W/m^2 """  ## not used
    return 2*h*nu_Hz**3/c2/expm1(h*nu_Hz/(kB*T_K))
    
def dBnu(nu_Hz, dT_K, T_K=Tcmb):
    """ Bnu(T+dT) - Bnu(T), black body, in W/m^2 """
    xT = h*nu_Hz/kB
    return 2*h*nu_Hz**3/c2 * (1/expm1(xT/(T_K+dT_K))-1/expm1(xT/T_K))

def cmb(nu_GHz, dT, T0_K=Tcmb, Omega=solid_angle):
    nu = nu_GHz * 1e9   ## Hz
    return 1e26*Omega*dBnu(nu, 1e-6*dT, T0_K)
    
def cmb_RJ(nu_GHz, dT_muK, Omega=solid_angle):
    nu = nu_GHz * 1e9   ## Hz
    return 1e-6*dT_muK*2.0e26 *kB*Omega*nu*nu/c2      ## is this right?  ## should use explicit equation above?
    
def synch(alpha, nu_GHz):
    return nu_GHz**alpha
    
def dust(tau_250, beta, T_dust, nu_GHz, Omega=solid_angle):
#     return tau_250*greybody(beta, T_dust, nu_GHz)*Omega
    nu = 1e9*nu_GHz  ## Hz
    x = h*nu/(kB*T_dust)
    S_dust = tau_250*2*(h*nu*nu*nu/c2)/expm1(x)*(nu/1.2e12)**beta
    return 1e26*S_dust*Omega   ## convert to Jy 

class M31model(submmModel2_normalized):
    nparam = 8
    fmtstring = "%.3f "*nparam
    paramBlocks =  range(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$\tau_{250}$", r"$\beta_{dust}$", r"$T_{dust}$", r"EM", r"$\Delta T_{CMB}$", r"$A_{synch}$", 
                r"$\alpha_{synch}$", r"$A_{AME}$"]
                
    # default prior params -- change with set_prior below
    EM_inv_sigma2 = 0.0   #(0.5)**(-2.0)
    EM_mean = 7.7
    alpha_synch_inv_sigma2 = 0.0 # (0.3)**(-2.0)
    alpha_synch_mean = -0.9


    def __init__(self, tau250, beta_dust, T_dust, EM, DT_CMB, A_synch, alpha_synch, A_ame):

        self.tau250 = tau250
        self.beta_dust = beta_dust
        self.T_dust = T_dust
        self.EM = EM
        self.DT_CMB = DT_CMB
        self.A_synch = A_synch
        self.alpha_synch = alpha_synch
        self.A_ame = A_ame
        

    @classmethod
    def set_prior(cls, EM_inv_sigma2, EM_mean, alpha_synch_inv_sigma2, alpha_synch_mean):
        if EM_inv_sigma2 > 0:
            print "setting EM prior to %f +- %f" % (EM_mean, EM_inv_sigma2**(-0.5))
            
        if alpha_synch_inv_sigma2 > 0:
            print "setting alpha_synch prior to %f +- %f" % (alpha_synch_mean, alpha_synch_inv_sigma2**(-0.5))
            
        cls.EM_inv_sigma2 = EM_inv_sigma2
        cls.EM_mean = EM_mean
        cls.alpha_synch_inv_sigma2 = alpha_synch_inv_sigma2
        cls.alpha_synch_mean = alpha_synch_mean

    @classmethod
    def prior(cls, tau250, beta_dust, T_dust, EM, DT_CMB, A_synch, alpha_synch, A_ame):

        if tau250<0 or EM<0 or A_synch<0 or A_ame<0:
            return 0
            
        if alpha_synch <-2.0 or alpha_synch>-0.5:
            return 0
            
        if DT_CMB<-200 or DT_CMB>200:  ### muK -- check units...
            return 0

        pri = 1
                    
        if cls.EM_inv_sigma2 > 0:
            pri = np.exp(-0.5*cls.EM_inv_sigma2*(EM-cls.EM_mean)**2)
            
        if cls.alpha_synch_inv_sigma2 > 0:
            pri *= np.exp(-0.5*cls.alpha_synch_inv_sigma2*(alpha_synch-cls.alpha_synch_mean)**2)
            
        return pri
            
            
    def at_nu(self, nu_GHz):
        return dust(self.tau250, self.beta_dust, self.T_dust, nu_GHz) + freefree(self.EM, nu_GHz) +\
               cmb(nu_GHz, self.DT_CMB) + self.A_synch * synch(self.alpha_synch, nu_GHz) + self.A_ame * AME(nu_GHz)
               
    def at(self, data):
        return dust(self.tau250, self.beta_dust, self.T_dust, data.freq) + freefree(self.EM, data.freq) +\
               cmb(data.freq, self.DT_CMB) + self.A_synch * synch(self.alpha_synch, data.freq) + self.A_ame * AME(data.freq)
               
    def all_at_nu(self, nu_GHz):
        return np.array([dust(self.tau250, self.beta_dust, self.T_dust, nu_GHz),
                         freefree(self.EM, nu_GHz),
                         cmb(nu_GHz, self.DT_CMB),
                         self.A_synch * synch(self.alpha_synch, nu_GHz),
                         self.A_ame * AME(nu_GHz)
                         ])
               
    __call__ = at    


### most of the following are generic and could be inherited
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

    def plot(self, data, logplot=True, wavelength=False, components=True, textfile=None):
        """plot the data and the model"""
        npt = 100
        if not logplot:
            f = linspace(min(data.freq), max(data.freq), npt)
        else:
            f = logspace(np.log10(min(data.freq)), np.log10(max(data.freq)), npt, 10)
        model_flux = self.at_nu(f)
        data.plot(fmt='o', wavelength=wavelength, logplot=logplot)
        if wavelength:
            f = speed_of_light/f
        plt.plot(f, model_flux)
        if components:            ## plot the individual components
            model_flux_comps = self.all_at_nu(f).transpose()
            plt.plot(f, model_flux_comps)
            ## legend???
            
        if textfile:
            header = "lambda" if wavelength else "nu"
            header += "\t" + "\t".join(['dust', 'ff', 'cmb', 'synch', 'AME'])
            np.savetxt(textfile, np.hstack((f.reshape(npt,1), model_flux_comps)), header=header)
            

        
    @classmethod
    def startfrom(cls, data=None, random=None):
        """
        generate a set of starting parameters for the model:
        tau250, beta_dust, T_dust, EM, DT_CMB, A_synch, alpha_synch, A_ame
        """
        cls.start_params = (1.0e-5, 2.0, 20., 8.0, 1.0, 10.0, -1.0, 50.0)  ## careful of units
        if random:
            ## these are all typically 1-2sigma 
            stds = (0.5e-5, 0.5, 1.5, 2.0, 2.0, 2.0, 0.3, 20)
            cls.start_params += np.random.randn(len(stds))*stds
            
        return cls.start_params
