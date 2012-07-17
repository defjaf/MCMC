
from __future__ import division

import math

import matplotlib.pyplot as plt
from numpy import array, exp, asarray, cos, sin, sqrt, float64, linspace, log, errstate, min, max, vectorize
from scipy import special,integrate
#import numexpr as ne

speed_of_light = 299792.458 ### micron GHz
c2 = speed_of_light**2
kB = 0.0  ### boltzmann
h = 0.0 ### Planck
Tcmb = 2.726   # CMB temperature (K)

import Proposal

from model import greybody

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

def freefree(nu, Te=8000.0, Omega=0.00382794):
    g_ff = np.log(4.955e-2/nu) + 1.5*np.log(Te)
    tau_ff = 3.014e-2 * Te**(-1.5) * g_ff
    Tff = Te*(1-np.exp(-tau_ff))    ### need to deal with small tau? see above...
    Sff = 2*kB*Tff * Omega * nu*nu/c2
    return Sff
    
    
def cmb(nu, Omega=0.00382794):
    Scmb = 2*kB*Omega*nu*nu/c2   ## is this right?  ## should use explicit equation above?
    
def synch(alpha, nu):
    return nu**alpha
    
def ame(nu):
    pass


class M31Model(object):
    nparam = 8
    fmtstring = "%.3f "*nparam
    paramBlocks =  range(nparam)    #### not right with different marginalization?
    nBlock = max(paramBlocks)+1
    texNames = [r"$\tau_{250}$", r"$\beta_{dust}$", r"$T_{dust}$", r"EM", r"$\Delta T_{CMB}$", r"$A_{synch}$", 
                r"\alpha_{synch}$", r"$A_{AME}$"]

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
    def prior(cls, A_dust, beta_dust, T_dust, EM, DT_CMB, A_synch, alpha_synch, A_ame):

        if A_dust<0 or EM<0 or A_synch<0 or A_ame<0:
            return 0
            
        if sync_amp <-2.0 or sync_amp>-0.5:
            return 0
            
        if DT_CMB<-200 or DT_CMB>200:  ### muK -- check units...
            return 0
            
    def at_nu(self, nu):
        return self.A_dust * greybody(self.b1, self.T1, nu) + EM * freefree(nu) +\
               DT_CMB * cmb(nu) + A_synch * synch(alpha_synch, nu) + A_ame * ame(nu)
               
    def at(self, data):
        return self.A_dust * greybody(self.b1, self.T1, data.freq) + EM * freefree(data.freq) + \
               DT_CMB * cmb(data.freq) + A_synch * synch(alpha_synch, data.freq) + A_ame * ame(data.freq)
               

    __call__ = at    