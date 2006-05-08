from __future__ import division
import Likelihood

# requires the Cl_model_class.__call__(self) return Cl[ell, 0:3]
# where 0:3 represents TT, TE, EE

class binnedClLikelihood(Likelihood.Likelihood):
    """
    represents a likelihood function for a binned set of C_l(XY)
    """

    def lnLike1(self, model_vals=None, data=None):
        return -data.calcLnLike(model_vals)  ## really returns -lnLike=chi2/2
    

    def lnNorm1():
        return 0.0
    

    def setModel(self, params):
        """set the model for the given parameters"""
        self.this_model = self.model(params)
        self.model_vals = self.this_model()

    def lnLike(self, params): 
        """
        calculate the sum of the lnLike for each of the datasets
        i.e., \sum_i ln p(data_i |params) up to a param-independent term,
        """
        self.setModel(params)
        return sum( [ self.lnLike1(self.model_vals, d) for d in self.data ] )
