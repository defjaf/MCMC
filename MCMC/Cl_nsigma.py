"""
module to calculate and plot the number of sigma for CMB data points
with respect to some model
"""




from .ClData import ClData
import string
import os.path
import math
import astropy.io.fits as pyfits
from pylab import *

from numpy import arange, array, float64, transpose, zeros, logical_and

filename = "data_list.txt"

def Cl_nsigma(WMAP=False, filename = filename, lmin=0,lmax=1500):

    names = {}
    namef = "data_names.txt"
    print("NAMES:")
    with open(namef, "r") as f:
        for line in f:
            name, fil = string.split(line)
            names[fil]=name
            print(fil, names[fil])
    
    
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
    ncol = len(colorstring)
    symstring = 'o+s*'
    figure(1)
    hold(True)
    for (iset, set) in enumerate(data):
        lnLike.append(set.calcLnLike(Cl))
            
        if (set.name == 'WMAP' and not WMAP) or (set.name != "WMAP" and WMAP):
            continue
            
        exptname.append(set.name)
        try:
            rng = set.Clpol_idx['TT']
        except KeyError:
            continue
            
        BP = array([set.getWinBandpower(j, Cl) for j in range(set.num_points)])

        diffs = set.getdelta(BP)/sqrt(set.var)
        diffs = diffs[rng[0]:rng[1]]

        ell1 =  set.ell[rng[0]:rng[1]]
        lidx = logical_and(ell1<lmax, ell1>lmin)
        ell.extend(ell1[lidx])

        nsig.extend(diffs[lidx])

        #expt.extend( (rng[1]-rng[0])*[iset] )
        expt.extend( len(lidx)*[iset] )

        try:
            lab=names[set.name]
        except:
            lab=set.name
        print(lab, set.name)
        plot(ell1[lidx], diffs[lidx], colorstring[iset%ncol]+symstring[iset//ncol], label=lab)

    # ### this only works for the current ordering (WMAP first)
    # if WMAP:
    #     legend([string.split(set.name, '_')[0] for set in data[0]])
    # else:
    #     legend([string.split(set.name, '_')[0] for set in data[1:]])
    legend()
        
    hold(False)

    
    figure(2)
    hold(True)
    hi = hist(nsig, bins=40)

    nsig=array(nsig, dtype=float64)
    m = nsig.mean()
    v = (nsig**2).mean() - m**2

    print('%f +- %f' % (m, sqrt(v)))
    gauss = max(hi[0])*exp(-0.5 * (hi[1]-m)**2/v)
    plot(hi[1], gauss)

    hold(False)
    return hi
    ### get a histogram of the number of sigma

if __name__ == "__main__":
    ret = Cl_nsigma()
