#### not used any more!!
### probably could make this more generic with 'repackager, unpackager' fns?
class Gaussian2DProposal(Proposal.GaussianProposal):
    """
    represents a Proposal Density for MCMC.
    for the 2D gaussian, use a symmetric gaussian proposal in all parameters?
    use angle = angle mod 2pi   (or really pi or even pi/2?)
        sigma = |sigma| ???
    for now, just use independent gaussians in each parameter?
    """

    def __init__(self, sigmas, start=None):
        """
        initialize proposal based on data or whatever
        """
        Proposal.GaussianProposal.__init__(self, unpackage(sigmas), start)
        
    def newParams(self):
        """
        just repackages the generic parameters for the specific model
        """
        return package(Proposal.GaussianProposal.newParams(self))
    
    def setStart(self, start):
        Proposal.GaussianProposal.setStart(self, unpackage(start))
        return
