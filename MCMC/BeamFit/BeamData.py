

from numpy import asarray, float64, log, where, nonzero, math
from numpy.random import uniform, normal

import pylab

from ..GaussianData import GaussianData

# AHJ: modify/extend/subclass for non-uniform noise
class BeamData(GaussianData):
    """
    represents measurements of a beam.
    assume data at positions given by arrays x[], y[]
    with noise variance sigma^2 (should work for array or scalar sigma)
    """

    def __init__(self, xarr, yarr, data, sigma, cts=None):
        GaussianData.__init__(self, data, sigma)
#        if self.n != yarr.size() or self.n != xarr.size():
#            raise GaussianData.DataSizeError(self.n, xarr.size(), yarr.size())
        if self.n != yarr.size or self.n != xarr.size:
            raise GaussianData.DataSizeError(self.n, xarr.size, yarr.size)
        
        self.x=asarray(xarr, float64)
        self.y=asarray(yarr, float64)

        if cts is not None:
            self.cts = cts

    def plot(self):
        """ plot the data """
        pylab.contour(self.x, self.y, self.d)
        
    def stats(self, sigcut=None, do_abs=False):
        """
        calculate (<x>, <y>, var x, var y), weighting by data (where +ve)

        if do_abs, use |data| instead
        
        if sigcut!=None, only use data with self.d > self.sigma
        """
        
        #if sigcut is None: sigcut = 0.0

        if do_abs:
            dd = abs(self.d)
        else:
            dd = self.d
            
        offset = dd.mean()    ## should really do this out of the center
        dd -= offset          ## IS THIS ACTUALLY SHIFTING THE self.d???

        if sigcut is not None:
            idx = nonzero(dd > sigcut*self.sig)
            dd = dd[idx]
            xx = self.x[idx]
            yy = self.y[idx]
        else:
            xx = self.x
            yy = self.y
            
        norm = dd.sum()
        Ex = (xx*dd).sum()/norm
        Ey = (yy*dd).sum()/norm
        Ex2 = (xx*xx*dd).sum()/norm
        Ey2 = (yy*yy*dd).sum()/norm
        Vx = Ex2 - Ex*Ex
        Vy = Ey2 - Ey*Ey

        return Ex, Ey, Vx, Vy


## used to be a subclass of BeamData, but easier this way?
def BeamSim(beam_model, n, sigma=1.0, amplitude=1.0,
            xarr=None, yarr=None, data=None, xrng=None, yrng=None):
    """
    simulated data for a beam.
    can generate positions xarr, yarr if needed, in ranges xrange, yrange
    requires callable object beam_model(xarr, yarr) to generate signal
    """
    if xrng is None: xrng=(-1,1)
    if yrng is None: yrng=(-1,1)
    if xarr is None and xrng is not None:
        xarr=uniform(xrng[0], xrng[1], n)
    if yarr is None and yrng is not None:
        yarr=uniform(yrng[0], yrng[1], n)

    if data is None and sigma is not None:
        data = amplitude*beam_model.atxy(xarr, yarr) + normal(0, sigma, n)

    return BeamData(xarr, yarr, data, sigma)


