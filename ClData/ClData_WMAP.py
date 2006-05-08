from __future__ import division

import os.path
import math
import copy

#from numarray import array, int32, float64, bool8, arange
from numpy import array, int32, float64, bool8, arange, rank

import ClData_CosmoMC
from ClData_CosmoMC import lmax
from WMAP_likelihood import WMAP_likelihood as WMAP

class ClData_WMAP(ClData_CosmoMC.ClData_CosmoMC):
    """
    represents WMAP-formatted data (ClTT, ClTE).
    
    NB this is only used for plotting purposes: the raw data is read
    internal to the WMAP_likelihood module
    """

    def __init__(self):
        self.readDataset()

    def readDataset(self):
        mapfile = ['wmap_binned_tt_powspec_yr1_v1p1.txt',
                   'map_binned_te_powspec_yr1_v1.txt']
        mapdir = 'cmb/misc-data/MAP/'
        homedir = os.path.expandvars('${HOME}/home')
        if not os.path.exists(homedir):
            homedir = os.path.expandvars('${HOME}')

        self.name = 'WMAP'

        self.Clpol_idx = {}

        self.ell = []
        self.Cl = []
        self.sig = []
        self.inc_pol = []

        self.has_pol = True
        self.ncls = 2
        self.beam_uncertain = False
        self.has_xfactors = False
        self.calib_uncertainty = 0.0

        i = 0


        ### get the TT, TE ranges
        for fname, ch_type in zip(mapfile, ['TT', 'TE']):

            if ch_type == 'TT':
                cols = [0, 3, 4]
            elif ch_type == 'TE':
                cols = [0, 1, 2]
            self.Clpol_idx[ch_type] = [i, i]
            WMAPname ='/'.join((homedir, mapdir, fname))
            for line in open(WMAPname):
                line = line.strip()
                if len(line)>0 and line[0] != '#':
                    splitline = line.split()
                    self.ell.append(int(splitline[cols[0]]))
                    self.Cl.append(float(splitline[cols[1]]))
                    self.sig.append(float(splitline[cols[2]]))
                    self.inc_pol.append(False)
                    i += 1
            self.Clpol_idx[ch_type][1] = i

        self.ell = array(self.ell, dtype=int32)
        self.Cl = array(self.Cl, dtype=float64)
        self.obs = self.Cl    ### nb NOT a copy
        self.sig = array(self.sig, dtype=float64)
        self.inc_pol = array(self.inc_pol, dtype=bool8)
        self.var = self.sig**2
        self.num_points = len(self.Cl)

    def getWinBandpower(self, ibin, Cl):
        """
        return the bandpower if observing Cl for WMAP.
        just assume flat in l(l+1)Cl.
        (use "real" data or binned?)
        """
        if ranl(Cl)==2:
            ClTT = Cl[0]
        elif rank(Cl)==1:
            ClTT=Cl
        else:
            raise MyErr('Cl shape wrong in getWinBandpower:', Cl.shape)

        for ch_type, lrange in self.Clpol_idx.iteritems():
            ell = arange(len(ClTT), dtype=float64)   # **** STOPPED HERE, NOT DONE
        win = ell*(ell+1)
        win /= sum(win)
        bandpower = dot(ClTT[:maxl],
                        self.window[i, 0, self.win_min[i]:maxl])
        if self.has_pol and rank(Cl)==2 and self.inc_pol[i]:
            bandpower += (Cl[1:num_cls, win_min[i], maxl+1]*
                          self.window[i, 1:num_cls, win_min[i], maxl+1]).sum()
            #for l in xrange(self.win_min[i], maxl):
            #    bandpower += dot(Cl[1:num_cls, l], self.window[i, 1:num_cls, l])            
        return bandpower


    def calcLnLike(self, cl, init_WMAP=[True]):
        """
        calculate -ln(likelihood) for WMAP data with respect to
        C_l = cl = [ClTT, ClTE, ClEE]  [not l(l+1)Cl/2pi]
        """

        ttFile = 'WMAP/likelihood/tt_diag.dat.gz'
        ttOffDiag = 'WMAP/likelihood/tt_offdiag.dat.gz'
        teFile = 'WMAP/likelihood/te_diag.dat.gz'
        teOffDiag = 'WMAP/likelihood/te_offdiag.dat.gz'

        twopi = 2.0*math.pi

        if init_WMAP[0]:
            if lmax<WMAP.WMAP_lmax_TT:
                raise MyErr('lmax not large enough for WMAP')
            stat = WMAP.WMAP_init(ttFile, ttOffDiag, teFile, teOffDiag)
            init_WMAP[0] = False  ## careful -- static var???!!!

        l = arange(WMAP.WMAP_lmax_TT+1)
        clTT = cl[0,:WMAP.WMAP_lmax_TT+1]*l*(l+1)/twopi

        lnLike = WMAP.WMAP_LnLike_TT(clTT)
        
        if self.has_pol:  ## allow turning off the polarization calculation
            clTE = cl[1,:WMAP.WMAP_lmax_TT+1]*l*(l+1)/twopi
            clEE = cl[2,:WMAP.WMAP_lmax_TT+1]*l*(l+1)/twopi
            lnLike += WMAP.WMAP_LnLike_TE(clTT, clTE, clEE)

        return -lnLike    ## expects negative for some reason...

        
