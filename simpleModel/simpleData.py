from __future__ import division
from numpy import asarray, float64, log
from numpy.random_array import uniform, normal   #check

import pylab

from GaussianData import GaussianData

# AHJ: modify/extend/subclass for non-uniform noise
class simpleData(GaussianData):
    """
    1d linear data with errors
    """
    
    def __init__(self, xarr, data, sigma):
        GaussianData.__init__(self, data, sigma)
        if self.n != xarr.size():
            raise GaussianData.DataSizeError(self.n, xarr.size())
        
        self.x=asarray(xarr, float64)

    def plot(self):
        """ plot the data """
        pylab.plot(self.x, self.d)

class simpleSim(simpleData):
    """
    simulated 1d linear data
    """
    
    def __init__(self, model, n, sigma=1.0, xarr=None, data=None, xrng=None):
        self.n = n
        if xrng is None: xrng=(-1,1)
        if xarr is None and xrng is not None:
            xarr=uniform(xrng[0], xrng[1], n)

        if data is None and sigma is not None:
            data = model.atx(xarr) + normal(0, sigma, n)

        BeamData.__init__(self, xarr, data, sigma)

    
