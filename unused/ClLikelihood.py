# subclasses for Gaussian, correlated gaussian,  offset-ln, WMAP, etc.?
# how to make the product of multiple datasets automatic?


import Likelihood

class GaussianLikelihood(Likelihood.Likelihood):
    """
    Represents the likelihood for CMB C_l data.
    """

    def lnLike(self, params):
        """return ln p(data|params) up to a param-independent term"""
        
        self.setModel(params)

        chi2 = ((self.model_vals-self.data.d)**2/self.data.sig2).sum()

        return -chi2/2
        
        
    def lnNorm(self):
        """
        return the parameter-independent part of the ln-likelihood.
        """
        norm = (self.data.n-1)*log(2*pi)
        norm2 = self.data.lnDetN + dNid
        
        return -0.5*(norm + norm2)
    
    __call__ = lnLike
