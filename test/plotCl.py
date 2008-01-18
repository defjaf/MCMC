""" read from datasets and plot """

from __future__ import division

from ClData import ClData
import string
import os.path
import math
import pyfits
from pylab import *

from numpy import arange, array, float64, transpose, zeros

colorstring='bgrcmyk'
symstring='.s+xD'   ## dots
symstring='os+xD'   ## circles

mapdir = 'cmb/misc-data/WMAP06/'
homedir = os.path.expandvars('${HOME}/home')
if not os.path.exists(homedir):
    homedir = os.path.expandvars('${HOME}')
mapdir = '/'.join( (homedir, mapdir) )

def plotCl(filename, plotLike=False, plotBP=False):
    
    mapf = '/'.join( (mapdir, 'models/cmb_04546021.fits') )
    mapd = pyfits.getdata(mapf)
    
    ell = arange(mapd.shape[0])
    norm = 1e6
    ClTT = array(norm**2*mapd.field(0))
    llCl = ClTT*ell*(ell+1)/(2*math.pi)
    Cl = array([ClTT, zeros(len(ClTT), float64), zeros(len(ClTT), float64)])
    
    data = ClData.getClData(filename)
    
    figure(1)
    hold(True)
    for iset, set in enumerate(data):
        rng = set.Clpol_idx['TT']
        col = colorstring[iset % len(colorstring)]
        sym = symstring[iset // len(colorstring)]
        plot(set.ell[rng[0]:rng[1]], set.Cl[rng[0]:rng[1]], col+sym)
    
    plot(ell, llCl)
    
    legend([string.split(set.name, '_')[0] for set in data])
    
    ## need to deal with errors separately like this...
    for iset, set in enumerate(data):
        rng = set.Clpol_idx['TT']
        col = colorstring[iset % len(colorstring)]
        sym = symstring[iset // len(colorstring)]
        errorbar(set.ell[rng[0]:rng[1]], set.Cl[rng[0]:rng[1]],
                 set.sig[rng[0]:rng[1]], ecolor=col, fmt=None)
        if plotBP and set.name != 'WMAP':
            bandpowers = array([
                set.getWinBandpower(i, Cl) for i in range(set.num_points)])
            plot(set.ell[rng[0]:rng[1]], bandpowers[rng[0]:rng[1]])
    
    xlabel(r'$\ell$')
    ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$')
    
    hold(False)
    
    if plotLike:
        
        lnlikelist = []
        lnlike = []
        ampl = arange(0.5, 1.5, 0.01)
        for amplitude in ampl:
            lnlikelist.append(set.calcLnLike(amplitude*Cl) for set in data)
            lnlike.append(sum(lnlikelist[-1]))
            print amplitude, lnlikelist[-1]
        
        lnlike=array(lnlike)
        lnlikelist=array(lnlikelist)
        print lnlike
        
        figure(2)
        hold(True)
        for set, like in zip(data, transpose(lnlikelist)):
            plot(ampl, like-min(like), label=string.split(set.name, '_')[0])
        
        plot(ampl, lnlike-min(lnlike), linewidth=3)
        
        hold(False)
        
        return lnlike
    
    else:
        return None
        
def main():
    return plotCl("data_list.txt", plotLike=True)
                    

if __name__ == "__main__":
    ret = main()
