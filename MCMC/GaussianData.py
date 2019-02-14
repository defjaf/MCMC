# need a different class for correlated gaussian data?



from numpy import asarray, float64, log, dot, transpose

class GaussianData(object):
    """
    more-or-less arbitrary gaussian data (d, sigma).  will likely need
    to subclass for any real application in order to hold ancillary
    information (C_l bins, beam maps positions, etc).
    """
    def __init__(self, data, sigma):
        self.d = asarray(data, float64)
        self.n = self.d.size
        self.sig = asarray(sigma, float64)    # convert even if scalar

        self.nsig = self.sig.size   # need asarry: scalar.size doesn't exist
        self.sig2 = self.sig*self.sig

        if self.nsig == self.n:       # only checks size, not shape...
            self.lnDetN = (log(self.sig2)).sum()
        elif self.nsig == 1:
            self.lnDetN = self.n*log(self.sig2)
        else:
            raise DataSizeError(self.n, self.nsig)
            
        self.data_chi2 = None
        

    def quadform(self, A=None, B=None):
        """
            calculate A^T N^-1 B or A^T N^-1 A or data^T N^-1 data
            works for vector A,B  [and probably matrix -- untested]
            NB. for this concrete class, N^-1 = diag(1/sig^2)
            
            Special case for if A==B==None, use A=B=self.data.
            (stored statically)
            
        """
        if A is None and B is None:
            if self.data_chi2 is None:
                # self.data_chi2 = (self.d*self.d/self.sig2).sum()
                self.data_chi2 = dot(self.d,self.d/self.sig2)
            return self.data_chi2
        elif B is None:
            return dot(A.transpose()/self.sig2,A)
            #return (A*A/self.sig2).sum()
        elif A is None:
            raise NotImplementedError ### Pick a better exception here!
        else:
            return dot(A.transpose()/self.sig2,B)
            # return (A*B/self.sig2).sum()

    def chi2(self, vec=None):
        """ 
            return the chi^2 (data-vec)^T N^-1 (data-vec)
            if vec is None or vec==0, just do data^T N^-1 data
        """ 
        
        if vec is None: 
            return self.quadform()
        else:
            return self.quadform(self.d-vec)
            
            
    
class DataSizeError(Exception):
    def __init__(self, *value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        