"""
module to calculate and plot the number of sigma for CMB data points
with respect to some model
"""

from __future__ import division
from ClData import ClData
import string
import os.path
import math
import pyfits
from pylab import *
#from numarray import arange, array, float64, Error, transpose, zeros
from scipy import arange, array, float64, transpose, zeros

def Cl_nsigma(WMAP=False):

    filename = "data_list.txt"
    
    mapdir = 'cmb/misc-data/MAP/'
    homedir = os.path.expandvars('${HOME}/home')
    if not os.path.exists(homedir):
        homedir = os.path.expandvars('${HOME}')
        
    mapdir = '/'.join( (homedir, mapdir) )
        
    data = ClData.getClData(filename)
    
    
    mapf = '/'.join( (mapdir, 'models/cmb_04546021.fits') )
    mapd = pyfits.getdata(mapf)
    
    ell = arange(mapd.shape[0])
    norm = 1e6
    ClTT = array(norm**2*mapd.field(0))
    Cl = array([ClTT, zeros(len(ClTT), float64), zeros(len(ClTT), float64)])

    lnLike = []
    nsig = []
    ell = []
    expt = []
    exptname = []
    colorstring='bgrcmyk'
   
    figure(1)
    hold(True)
    for (iset, set) in enumerate(data):
        lnLike.append(set.calcLnLike(Cl))
            
        if (set.name == 'WMAP' and not WMAP) or (set.name != "WMAP" and WMAP):
            continue
            
        exptname.append(set.name)
        rng = set.Clpol_idx['TT']
        BP = array([set.getWinBandpower(j, Cl) for j in range(set.num_points)])

        diffs = set.getdelta(BP)/sqrt(set.var)
        diffs = diffs[rng[0]:rng[1]]
        nsig.extend(diffs)

        ell1 =  set.ell[rng[0]:rng[1]]
        ell.extend(ell1)

        expt.extend( (rng[1]-rng[0])*[iset] )

        plot(ell1, diffs, colorstring[iset]+'o', label=set.name)

    ### this only works for the current ordering (WMAP first)
    if WMAP:
        legend([set.name for set in data[0]])
    else:
        legend([set.name for set in data[1:]])
        
    hold(False)

    
    figure(2)
    hold(True)
    hi = hist(nsig, bins=40)

    nsig=array(nsig, dtype=float64)
    m = nsig.mean()
    v = (nsig**2).mean() - m**2

    print '%f +- %f' % (m, sqrt(v))
    gauss = max(hi[0])*exp(-0.5 * (hi[1]-m)**2/v)
    plot(hi[1], gauss)

    hold(False)
    return hi
    ### get a histogram of the number of sigma

if __name__ == "__main__":
    ret = Cl_nsigma()
