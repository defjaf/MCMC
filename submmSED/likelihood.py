from __future__ import division
import Likelihood

from numpy import linalg, dot


class SEDLikelihood1(Likelihood.Likelihood):
    """
    submm SED flux likelihood, marginalized over individual amplitudes.
    (This is essentially the same as the beamFit Likelihood.Likelihood, but boosted to matrices.)
    """
    ## this doesn't need to be a class method, although can just use self.model_vals!
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
        
        ## solve for (FNiF) z = (FNid) for z
        z = linalg.solve(FNiF, FNid)
        
        return 0.5 * (dot(FNid.transpose(), z) - log(detFNiF))


    def lnLike(self, params): 
        """
        calculate the sum of the lnLike for each of the datasets
        i.e., \sum_i ln p(data_i |params) up to a param-independent term
        
        assumes no correlations between datasets (which means between objects)
        """
        self.setModel(params)
        return sum( [ self.lnLike1(self.model_vals, d) for d in self.data ] )



class SEDLikelihood1(Likelihood.Likelihood):
    """submm SED flux likelihood, marginalized over overall amplitude.
       nb. need to be able to 
           a) have parameter r12 (ratio of component amplitudes) be object-specific
           b) fix those per-object parameters to 0 for single-component model
           
           but nb. the beamfit Likelihood.Likelihood.lnLike1 is exactly the correct form, except for 
           the object-specific part. So really the only change is in the full lnLike which needs to 
           set individual r12
    """


    def lnLike(self, params): ### NOT DONE
        """
        calculate the sum of the lnLike for each of the datasets
        i.e., \sum_i ln p(data_i |params) up to a param-independent term,
        """
        self.setModel(params)
        return sum( [ self.lnLike1(self.model_vals, d) for d in self.data ] )
