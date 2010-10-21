#!/usr/bin/env python
# encoding: utf-8
"""
blackbody.py

Created by Andrew H. Jaffe on 2010-09-22.
Copyright (c) 2010 Imperial College London. All rights reserved.

Module to calculate blackbody and greybody fluxes. 

Refactored out for conversion to cython
"""

import numpy as np
cimport numpy as np

from scipy import integrate, special
from scipy cimport integrate, special

cdef float h_over_k = 0.04799237 ###  K/Ghz
cdef float prefac = 1.0e-2 #### FIXME: find a physical definition to go here
cdef float nu0 = 100.0

expmax = 750.0


## note that nu is at least possibly a numpy array of floats

## the grey body function is (2h/c^2) (v/v0)^b v^3/(e^x-1)    w/x=hv/kT
## the integral of the grey body function is 
###(2h/c^2) (kT/h)^4 (kT/(hv0))^b Gamma[4 + b] Zeta[4 + b]
#### Gamma[4]=6, Zeta[4]=pi^4/90

## TODO: refactor blackbody (and greybody) for use as cython?
def blackbody(T, np.ndarray[np.float64_t, ndim=1] nu):
    """return the blackbody flux at frequency nu for temperature T [CHECK UNITS]"""
    
    x = h_over_k*nu/T
    ret = np.zeros_like(nu)
    idx = x<=expmax
    ret[idx] = prefac*nu[idx]**3/np.exp(x[idx]-1)
#     if x>expmax:
#         return 0.0
#     else:
# #    with np.errstate(over='ignore'):
#         return prefac*nu**3/(np.exp(x)-1)
# 
def greybody(beta, T, np.ndarray[np.float64_t, ndim=1] nu):
    """return the greybody flux at frequency nu for temperature T [CHECK UNITS]"""

    x = h_over_k*nu/T
    ret = np.zeros_like(nu)
    idx = x<=expmax
    ret[idx] = prefac*nu0**(-beta) * nu[idx]**(3+beta)/np.exp(x[idx]-1)
  #   x = h_over_k*nu/T
  #    if x>expmax:
  #        return 0.0
  #    else:
  # #   with np.errstate(over='ignore'):
  #      #  return ne.evaluate("prefac*nu**3/(exp(x)-1)")
  #      return prefac*nu**(3+beta)/(np.exp(x)-1)
  
  def totalflux(beta, T, nu1=None, nu2=None):
    """
    calculate the total flux of a grey body (with prefactors defined as above) over (nu1,nu2)
    """

    if nu1 is None and nu2 is None:
        ## return analytic expression
        return prefac * nu0**(-beta) * (T/h_over_k)**(4+beta) * \
               special.gamma(4+beta)*special.zeta(4+beta,1)
    else:
        ## do numeric integral
        return integrate.quad(lambda nu: greybody(beta, T, nu), nu1, nu2)
        
    