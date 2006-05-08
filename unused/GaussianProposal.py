
### probably unused; see genericGaussianProposal below
class GaussianProposal(object):
    """
    represents an uncorrelated multivariate Proposal Density for MCMC
    """

    def __init__(self, sigmas, start=None):
        """
        initialize proposal based on data or whatever
        """
        self.sigmas = sigmas
        self.n = sigmas.len()
        if start is not None: self.setStart(start)
        return

    def setStart(self, params):
        self.newParams = params
        return

    ## do this so we can package things differently in subclasses?
    def getSigmas(self):
        return self.sigmas

    def getCurrentParams(self):
        return self.newParams

    def getNewParams(self, prev):
        """
        get a new set of parameters based on the previous
        """
        self.prevParams = self.newParams.copy
        self.newParams += normal(0,1, self.n)*self.getSigmas()
        return self.newParams

    ## don't really need this for this symmetric case!
    ## ignore? or take advantage of fact that it's f(|x-y|)
    def density(self, params1, params2):
        """
        get the value of the density at a set of params
        (only needed for MCMC is proposal is asymmetric in (old, new)
        """
        pass

    def densityAtNew(self):
        """
        get the value of the density at the newly-generated parameters
        """
        return self.density(self.newParams,self.prevParams)

    def densityRatio(self):
        """
        get the density ratio q(old, new)/q(new, old) needed for
        Metropolis-Hastings.
        Just 1 if symmetric: q(old, new)=q(new, old)
        """
        return 1
