import GaussianData

class ClData0(GaussianData.GaussianData):
    """
    C_l data:
       array of class ClDataPoint for each bin
          this may not make sense for numarray-type calculations;
          don't want to convert each of these to an array each time
          (maybe just important for <C>, variance, x?)
       optional correlation matrix
       (if x_bin not present, assume gaussian??)
    """
    def __init__(self, lmin, lmax, C_bin, sigma_bin,
                 x_bin=None, windows=None):
        GaussianData.__init__(self, C_bin, sigma_bin)

    def readCosmoMCData(self, datafile):
        pass


### nb. want to be able to access each of these as an array for speed
#       so probably just package this up as method of ClData (or DataSet?)?
class ClDataPoint(object):
    """
    A single CMB data point:
       lmin, lmax
       window function from lmin -> lmax
           store as an array win[i] = window(lmin+i)?
       <C>, err_min, err_pls, sigma, var
       BJK x factor
       beam error
    """
    def __init__(lmin=None, lmax=None, var=None):
        if lmin is None:
            lmin = 0
        if lmax is None:
            lmax = 0
        self.range = (lmin, lmax)
        self.var = var

    def win(self, ell):
        return self.winarr[ell-lmin]

