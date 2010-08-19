from __future__ import division
import Likelihood

from numpy import linalg, dot, log, empty, log10

from Likelihood import ZeroPosterior

dolog10 = False  # "compile-time" flag for the derived parameter

class SEDLikelihood2(Likelihood.Likelihood):
    """
    submm SED flux likelihood, marginalized over individual amplitudes.
    (This is essentially the same as the beamFit Likelihood.Likelihood, but boosted to matrices.)
    
    Current version assumes the same parameters for each object 
    (so to do a single object just use a data-list of length 1)
    
    
    """
    ## this doesn't need to be a class method, although can just use self.model_vals!
    
    ## DONE appear to be problems here -- returning negative likelihood for reasonable values.
    ##      needed better normalization in blackbody function -- FIXME: make physical
    ##      is there any way to deal with a purposely-zero eigenvalue in the determinant
    ##        (e.g., when one amplitude is fixed to zero?)
    
    def lnLike1(self, model_vals=None, data=None):
        """
        return ln p(data|params) up to a param-independent term,
        - for a single dataset
        - assuming params have already been set by self.setModel
        
        multivariate extension of Likelihood.Likelihood formulae:
                if F = BA + d    <dd>=N
        marginalized over A, 
        -2lnP = ln | N (BB) | + (FF) - (FB)(BB)^-1(BF)
           w/ (BB) = B^T N^-1 B etc.
        """
        ## nb. FNid = F^T N^-1 d, etc.

        FNiF = data.quadform(model_vals)     
        FNid = data.quadform(model_vals, data.d)
        
        detFNiF = linalg.det(FNiF)
        if detFNiF<=0:
            raise ZeroPosterior(FNiF)
            
        #### want (FNid)^T (FNiF)^-1 (FNid)
        
        ## solve for (FNiF) z = (FNid) for z = (FNiF)^{-1} (FNid)
        ## NB. z is the ML value of the marginalized amplitudes!!!
        self.MLamplitude = linalg.solve(FNiF, FNid)

        # ## require positive amplitudes
        # if self.MLamplitude<0:
        #     raise ZeroPosterior(self.MLamplitude)
                        
        return 0.5 * (dot(FNid.transpose(), self.MLamplitude) - log(detFNiF))
        
    nDerived = 2  ## actually it depends on the dimension of, e.g., model_vals

    if dolog10:
        derivedTexNames = [r"$\ln A_1$", r"$\ln A_2$"]
        def getDerived(self, *params):
            """ calculate a list of any derived parameters """  
            return log10(self.MLamplitude)
    else:
        derivedTexNames = [r"$A_1$", r"$A_2$"]
        def getDerived(self, *params):
            """ calculate a list of any derived parameters """  
            return self.MLamplitude
        
        
        

class SEDLikelihood1(Likelihood.Likelihood):
    """
       single-component grey-body submm SED flux likelihood, marginalized over overall amplitude.
    """

    nDerived = 1    
    if dolog10:
        derivedTexNames = [r"$\ln A$"]
        def getDerived(self, *params):
            """ calculate a list of any derived parameters """  
            return [log10(self.MLamplitude)]
    else:
        derivedTexNames = [r"$A$"]
        def getDerived(self, *params):
            """ calculate a list of any derived parameters """  
            return [self.MLamplitude]


        
        
class SEDLikelihood_normalized(Likelihood.Likelihood):
    """
    submm SED flux likelihood, including individual amplitudes.
    (This is essentially a Gaussian likelihood)

    Current version assumes the same parameters for each object 
    (so to do a single object just use a data-list of length 1)
    
    Should work for both 1- and 2-parameter models

    """
    nDerived = 0
    
    def lnLike1(self, model_vals=None, data=None):
        """
        return ln p(data|params) up to a param-independent term,
        - for a single dataset
        - assuming params have already been set by self.setModel

        Gaussian likelihood
            if d_i = A B(p_i) + e_i    <ee>=N
            
        -2 ln P = (d - AB)^T N^{-1} (d - AB) + ln det N
                = chi^2 + ln det N
        """
        ## nb. dNiB = d^T N^-1 B, etc.

        return -0.5 * data.chi2(model_vals)


    def lnNorm1(self, data=None):
        """
        return the parameter-independent part of the ln-likelihood
        - for a single dataset
        """

        return -0.5 * ( data.n * log(2*pi) + data.lnDetN)
        


