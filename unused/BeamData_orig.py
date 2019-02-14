
from numarray import asarray, float64
from numarray.random_array import uniform, normal

# AHJ: modify/extend/subclass for non-uniform noise

class BeamData(object):
    """
    represents measurements of a beam.
    assume unbinned data at positions given by arrays x[], y[]
    with noise variance sigma^2 (should work for array or scalar sigma)
    """

    def __init__(self, xarr, yarr, data, sigma):
        self.n = xarr.size()
        self.nsig = len(asarray(sigma))   # need asarry since len(scalar) doesn't exist
        if (self.n != yarr.size() or self.n != data.size() or 
            self.nsig != 1 and self.nsig != self.n):
            raise DataSizeError(self.f, self.nsig)
        
        self.x=asarray(xarr, float64)
        self.y=asarray(yarr, float64)
        self.d=asarray(data, float64)
        self.sig=sigma
        self.sig2=sigma*sigma

        if self.nsig == n:
            self.lnDetN = (log(self.sig2)).sum()
        else:
            self.lnDetN = self.n*log(self.sig2)

        

class BeamSim(BeamData):
    """
    simulated data for a beam.
    can generate positions xarr, yarr if needed, in ranges xrange, yrange
    requires callable object beam_model(xarr, yarr) to generate signal
    """
    
    def __init__(self, beam_model, n, sigma=1.0, amplitude=1.0,
                 xarr=None, yarr=None, data=None, xrng=None, yrng=None):
        self.n = n
        if xrng is None: xrng=(-1,1)
        if yrng is None: yrng=(-1,1)
        if xarr is None and xrng is not None:
            xarr=uniform(xrng[0], xrng[1], n)
        if yarr is None and yrng is not None:
            yarr=uniform(yrng[0], yrng[1], n)

        if data is None and sigma is not None:
            data = amplitude*beam_model.at2(xarr, yarr) + normal(0, sigma, n)

        BeamData.__init__(self, xarr, yarr, data, sigma)


class DataSizeError(Exception):
    def __init__(self, *value):
        self.value = value
    def __str__(self):
        return repr(self.value)

