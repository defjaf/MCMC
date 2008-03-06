from __future__ import division

from numpy.random import normal, seed
from numpy import dot, where, zeros, asarray, float64
import numpy.linalg as la
from copy import copy

#debugf = file("debug.out", "w")


#TODO: modify to allow rotated parameters based on some correlation matrix
#      add flag self.rotate to __init__? or subclass?

class GenericGaussianProposal(object):
    """
    multivariate Proposal Density for MCMC

    if uncorrelated, sigmas=std devn
    if correlated,
         sigmas=sqrt(diagonal elements of matrix) and
         normalized_matrix=M(i,j)/[sig(i)*sig(j)]
       or, can just leave sigmas=ones(n) and normalized_matrix=M

       can specify M, M^{1/2}, depending on values of 
       matrix, isSqrt
                 
    Requires (a model class with?) package and unpackage methods which convert
    to and from the model's 'native' parameter format. If None, default to
    the identity
    
    allow 'fixed' parameters via sigma=0 proposal density
    
    allow updating only sets of parameters, defined by 
       self.block = [list of integers 0...nblocks-1]

    """

    def __init__(self, sigmas=None, package=None, unpackage=None,
                 matrix=None, isSqrt=False, rotateParams=False):
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
        
        seed(seed=None)

        if sigmas is not None:
            self.setSigmas(sigmas)

        ### or do this with classes???
        if matrix is not None:
            if not isSqrt:
                self.setNormalizedMatrix(matrix)
            else:
                self.setNormalizedMatrixFromSqrt(matrix)
        else:
            self.sqrtMatrix = None
            
        self.rotateParams=rotateParams
        

    def setNormalizedMatrix(self, matrix):
        """ set the Normalized correlation matrix

        """
        ## needs to excise the rows/cols with sigmas==0 before the cholesky
        #### COULD DO THIS WITH A MASK!!!???

        np = len(self.subset)
        if self.n == np:
            mat2 = matrix
        else:
            mat2 = [matrix[i,j] for i in self.subset for j in self.subset]
            mat2 = asarray(mat2).reshape((np, np))


        sqrtMatrix=la.cholesky(mat2)
        
        self.sqrtMatrix = zeros((self.n, self.n), dtype=float64)
        for i in range(np):   ### there's certainly a fancy-indexing version of this which works...
            for j in range(np):
                self.sqrtMatrix[self.subset[i],self.subset[j]] = sqrtMatrix[i,j]
                
        # DEBUG
#        print "set sqrtMatrix=\n", self.sqrtMatrix

    def setNormalizedMatrixFromSqrt(self, matrix):
        """ set the Normalized correlation matrix from the sqrt
        """
        self.sqrtMatrix=matrix
        # DEBUG
#        print "set sqrtMatrix=", self.sqrtMatrix
          
          
    def setSigmas(self, sigmas):
        self.sigmas = asarray(self.unpackage(sigmas))
        self.n = len(self.sigmas)
        self.subset = where(self.sigmas>0)[0]    # nb. where() returns a tuple
        
        ## DEBUG
#        print "set self.sigmas=", self.sigmas
#        print "    self.n=", self.n
#        print "    self.subset=", self.subset

    def getSigmas(self):
        return self.package(self.sigmas)

    def getNewParams(self, prevParams, block=None):
        """
        get a new set of parameters based on the previous
        - allow 'fixed' parameters via sigma=0 proposal density
        - only change parameters where the index is in block
        
        if rotateParams is set, sample from orthogonal parameters (but output "physical" parameters)
        
        if block is set, sample the specific parameter[s] in the block
        (nb. this is true even when rotateParams is set)
        
        """

        ### nb. copy crucial here otherwise you just change the prevParams
        ### when unpackage doesn't create a new object (e.g., the identity)!

        self.newParams = copy(self.unpackage(prevParams))
        
        ### need to check shape of block
        
        if self.rotateParams:
            assert block is not None and self.sqrtMatrix is not None
            offset = sum(self.sqrtMatrix[:,j]*normal(0,1) for j in block)
            self.newParams += offset*self.sigmas
        else:
            if self.sqrtMatrix is None:
                offset = normal(0,1, self.n)*self.sigmas
            else:
                offset = dot(self.sqrtMatrix,
                             normal(0,1, self.n))*self.sigmas                           

            if block is None:
                self.newParams += offset
            else:
                self.newParams[block] += offset[block]
            
        ## DEBUG
#        print >> debugf, block, offset

        return self.package(self.newParams)

    ## don't really need this for this symmetric case!
    ## ignore? or take advantage of fact that it's f(|x-y|)
    def value(self, params1, params2):
        """
        get the value of the proposal density at a set of params
        (only needed for MCMC is proposal is asymmetric in (old, new)
        """
        raise NotImplementedError


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
        Just ln1=0 if symmetric: q(old, new)=q(new, old)
        """
        return 0
   
