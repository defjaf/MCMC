from __future__ import division
import Likelihood

from numpy import linalg, dot, log, empty

from Likelihood import ZeroPosterior


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
    derivedTexNames = [r"$\ln A_1$", r"$\ln A_2$"]
    def getDerived(self, model_vals=None, data=None):
        """ calculate a list of any derived parameters """  
        return log(self.MLamplitude)
        
        

class SEDLikelihood1(Likelihood.Likelihood):
    """
       single-component grey-body submm SED flux likelihood, marginalized over overall amplitude.
    """

    nDerived = 1    
    derivedTexNames = [r"$\ln A$"]
    def getDerived(self, model_vals=None, data=None):
        """ calculate a list of any derived parameters """  
        return [log(self.MLamplitude)]
        


        
        
    
    
    