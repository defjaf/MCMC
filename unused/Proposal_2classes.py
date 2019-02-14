
from numarray.random_array import normal
from numarray import matrixmultiply
import numarray.linear_algebra as la

class UncorrelatedGaussianProposal(object):
    """
    multivariate Proposal Density for MCMC

    if uncorrelated,  sigmas=std devn

                 
    Requires (a model class with?) package and unpackage methods which convert
    to and from the model's 'native' parameter format. If None, default to
    the identity
    """

    def __init__(self, sigmas=None, package=None, unpackage=None,
                 matrix=None, isSqrt=False):
        """
        initialize proposal based on appropriate variances

        see class doc for how to specify correlations
        
        for models in which the parameters aren't a flat sequence need
        to supply 'package' and 'unpackage' operations to convert
        to/from the native format
        """

        if package is None:
            package = (lambda x:x)   #identity operator
        if unpackage is None:
            unpackage = (lambda x:x)

        self.package=package      #nb these are essentially staticmethods
        self.unpackage=unpackage

        if sigmas is not None:
            self.setSigmas(sigmas)


    def setSigmas(self, sigmas):
        self.sigmas = self.unpackage(sigmas)
        self.n = len(self.sigmas)        

    def getSigmas(self):
        return self.package(self.sigmas)

    def getNewParams(self, prevParams):
        """
        get a new set of parameters based on the previous
        """
        
        self.newParams = (self.unpackage(prevParams) + 
                          normal(0,1, self.n)*self.sigmas)
        
        return self.package(self.newParams)

    ## don't really need this for this symmetric case!
    ## ignore? or take advantage of fact that it's f(|x-y|)
    def value(self, params1, params2):
        """
        get the value of the proposal density at a set of params
        (only needed for MCMC is proposal is asymmetric in (old, new)
        """
        raise "Should never be here"


    #### get rid of self.new*, self.prev*???
    def densityAtNew(self):
        """
        get the value of the density at the newly-generated parameters
        """
        return self.value(self.unpackage(self.newParams),
                          self.unpackage(self.prevParams))

    def densityRatio(self, old, new):
        """
        get the density ratio q(old, new)/q(new, old) needed for
        Metropolis-Hastings.
        Just 1 if symmetric: q(old, new)=q(new, old)
        """
        return 1
   
    def lndensityRatio(self, old, new):
        """
        get the ln (density ratio) =ln q(old, new)/q(new, old) needed for
        Metropolis-Hastings.
        Just 1 if symmetric: q(old, new)=q(new, old)
        """
        return 0
   

class CorrelatedGaussianProposal(UncorrelatedGaussianProposal):
    """
    Correlated Gaussian proposal density.

    like uncorrelated proposal density, except
    
    sigmas=sqrt(diagonal elements of matrix) and
    normalized_matrix=M(i,j)/[sig(i)*sig(j)]
    or, can just leave sigmas=ones(n) and normalized_matrix=M
    
    can specify M, M^{1/2}, depending on values of 
    matrix, isSqrt
    """
    
    def __init__(self, sigmas=None, package=None, unpackage=None,
                 matrix=None, isSqrt=False):
        """
        initialize proposal based on appropriate covariances

        see class doc for how to specify correlations
        
        """
        
        UncorrelatedGaussianProposal.__init__(self, sigmas=sigmas,
                                              package=package, unpackage=unpackage)

        if matrix is not None:
            if not isSqrt:
                self.setNormalizedMatrix(matrix)
            else:
                self.setNormalizedMatrixFromSqrt(matrix)
        else:
            self.sqrtMatrix = None

    def setNormalizedMatrix(self, matrix):
        """ set the Normalized correlation matrix 
        """
        self.sqrtMatrix=la.cholesky_decomposition(matrix)

    def setNormalizedMatrixFromSqrt(self, matrix):
        """ set the Normalized correlation matrix from the sqrt
        """
        self.sqrtMatrix=matrix


    def getNewParams(self, prevParams):
        """
        get a new set of parameters based on the previous
        """

        self.newParams = self.unpackage(prevParams)
        self.newParams += matrixmultiply(self.sqrtMatrix,
                                         normal(0,1, self.n)*self.sigmas) 
        
        return self.package(self.newParams)    
