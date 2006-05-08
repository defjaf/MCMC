from __future__ import division
import Likelihood

class NormalizedBeamLikelihood(Likelihood.Likelihood):
    """
    likelihood function for a normalized beam
    """

    def lnLike1(self, model_vals, data):
        """
            return the log of the likelihood for the normalized beam model
            (where 'normalized' here means that we do care about the 
            normalization; we don't marginalize over it)
        """
        ## just return the chi^2
        return -0.5*data.chi2(model_vals)

    def lnNorm1(self, data=None):
        return -0.5*data.lnDetN

