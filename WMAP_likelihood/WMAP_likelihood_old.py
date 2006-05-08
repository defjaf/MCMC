"""
 WMAP likelihood code
 Written by Licia Verde and Hiranya Peiris, Princeton University,
      December 2002.

 F90 version by Antony Lewis Feb 03
 Python version by Andrew Jaffe Aug 05
"""

from __future__ import division
from numarray import array, Float64, sqrt, log, sum, zeros, arange, where


### module variables
WMAP_lmax_TT = 900
WMAP_lmax_TE = 450
WMAP_lmax_TE_file = 512 

WMAP_lmin_TT = 2
WMAP_lmin_TE = 2     #Can change l_min here to remove e.g. quadrupole 

num_diag_TT = (WMAP_lmax_TT-1)*(WMAP_lmax_TT-2)//2
num_diag_TE = (WMAP_lmax_TE-1)*(WMAP_lmax_TE-2)//2

# TT
cl_data = array(shape=WMAP_lmax_TT+1, type=Float64)
neff = array(shape=WMAP_lmax_TT+1, type=Float64)
fskyeff = array(shape=WMAP_lmax_TT+1, type=Float64)
r_off_diag = array(shape=num_diag_TT, type=Float64)
off_diag = array(shape=num_diag_TT, type=Float64)

#TE data
te_data = array(shape=WMAP_lmax_TE+1, type=Float64)
ntt = array(shape=WMAP_lmax_TE+1, type=Float64)
nee = array(shape=WMAP_lmax_TE+1, type=Float64)
te_off_diag = zeros(shape=num_diag_TE, type=Float64)
te_tt = array(shape=WMAP_lmax_TE+1, type=Float64)

def WMAP_init_TT (clFile, offDiag):

    fp = open(clFile)
    for l in xrange(2, WMAP_lmax_TT+1):
        line = fp.readline().split()
        try:
            if len(line):
                cl_data[l], neff[l], fskyeff[l] = [float(e) for e in line[1:4]]
        except:
            print 'wmap_init_tt:', l, line

    ix = 0
    fp = open(offDiag)
    for l in xrange(2, WMAP_lmax_TT+1):
        for ll in xrange(l+1, WMAP_lmax_TT+1):
            line = fp.readline().split()
            if len(line):
                i, j = int(line[0]), int(line[1])
                off_diag[ix], r_off_diag[ix] = float(line[2]), float(line[3])
                if l >= WMAP_lmin_TT:
                    ix+=1
                if l!=i or ll!=j:
                    raise MyErr('Error reading TT off diag (%d %d), (%d %d)'%( l, i, ll, j))

    return 0  ## should return an error code


def WMAP_init_TE (clFile, offDiag):

    fp = open(clFile)
    for l in xrange(2, WMAP_lmax_TE+1):
        line = fp.readline().split()
        if len(line):
            te_data[l], te_tt[l], ntt[l], nee[l] =[float(e) for e in line[1:5]]

    ix = 0
    fp = open(offDiag)
    for l in xrange(2, WMAP_lmax_TE+1):
        for ll in xrange(l+1, WMAP_lmax_TE_file+1):
            line = fp.readline().split()
            if len(line):
                i, j = int(line[0]), int(line[1])
                tmp = float(line[2])
                if l!=i or ll!=j:
                    raise MyErr('Error reading TE off diag')
                if l >= WMAP_lmin_TE and ll<WMAP_lmax_TE:
                    te_off_diag[ix] = tmp
                    ix+=1

    return 0  ## should return an error code

def WMAP_init(TclFile, ToffDiag, TEclFile, TEoffDiag):
    stat = WMAP_init_TT(TclFile, ToffDiag)
    if stat != 0:
        return stat
    return WMAP_init_TE(TEclFile, TEoffDiag)

def WMAP_LnLike_TT(clth):
    """
    WMAP TT likelihood: modified AHJ Aug 05 to return
    log(1e-10) if (cl+neff)<0
    """
    
    lmax = min(WMAP_lmax_TT+1, len(clth))  ## really l_max+1
    l = arange(WMAP_lmin_TT, lmax)

    Fdiag = zeros(lmax, type=Float64)
    Fdiagsqrt = zeros(lmax, type=Float64)
    z = zeros(lmax, type=Float64)
    zbar = zeros(lmax, type=Float64)

    ## nb global objects are indexed by l=lmin:lmax+1,
    ## local are indexed: 0:lmax-lmin+1
    Fdiag[l] = 2.0*(clth[l]+neff[l])**2 / ((2.0*l+1.0)*fskyeff[l]**2)
    Fdiagsqrt[l] = 1.0/sqrt(Fdiag[l])

    zth = clth[l] + neff[l]
    zth[where(zth<=0.0)] = 1.0e-10   ### AHJ

    z[l] = log(cl_data[l]+neff[l])
    zbar[l] = log(zth)

    Fisher = 1/Fdiag[l]
    off_log_curve = (clth[l]+neff[l])**2 * Fisher
    dchisq = 2.0/3.0 * (z[l]-zbar[l])**2 *off_log_curve + \
             1.0/3.0 * (clth[l]-cl_data[l])**2 * Fisher

    chisq = sum(dchisq)

    if chisq < WMAP_lmax_TT*2:
        #Only get off-diagonal terms if not a really bad fit, otherwise they
        #will be wildly wrong
        offchisq = 0
        ix = 0
        for l in xrange(WMAP_lmin_TT, lmax):
            for ll in xrange(l+1, lmax):
                Fisher = r_off_diag[ix]*Fdiagsqrt[l]*Fdiagsqrt[ll] \
                         + off_diag[ix]/(Fdiag[l]*Fdiag[ll])
                off_log_curv =(clth[l] + neff[l])*Fisher*(clth[ll] + neff[ll])
                dchisq = 2./3.*(z[l]-zbar[l])*off_log_curv*(z[ll]-zbar[ll]) + \
                         1./3.*(clth[l]-cl_data[l])*Fisher*(clth[ll]-cl_data[ll])
                offchisq += dchisq
                ix += 1

        chisq += 2*offchisq

    return -chisq/2.0
    
                
def WMAP_LnLike_TE(cltt, clte, clee):
    fsky = 0.85
    chisq = 0.0

    lmax = min(WMAP_lmax_TE+1, len(cltt), len(clte), len(clee))
    Fdiagsqrt = zeros(lmax, Float64)



    l = arange(WMAP_lmin_TE, lmax, type=Float64)
    ztt = cltt[l] + ntt[l]
    ztt[where(ztt<=0.0)] = 1.0e-10   ### AHJ
    zee = clee[l] + nee[l]
    zee[where(zee<=0.0)] = 1.0e-10   ### AHJ

    Fdiag = (ztt*zee + clte[l]*clte[l])/((2.0*l+1.0)*fsky**2/1.14)
    Fdiagsqrt[l] = 1/sqrt(Fdiag)

    delta_chisq = (clte[l] - te_data[l])**2/Fdiag

    chisq = sum(delta_chisq)

    offchisq = 0
    ix = 0
    for il in xrange(WMAP_lmin_TE, lmax):
        for ill in xrange(il+1, lmax):
            Fisher = te_off_diag[ix]*Fdiagsqrt[il]*Fdiagsqrt[ill]
            delta_chisq = (clte[il] - te_data[il])*Fisher* \
                          (clte[ill]-te_data[ill])
            offchisq +=delta_chisq
            ix += 1

    chisq += offchisq*2.0

    return -chisq/2.0


class MyErr(Exception):
    def __init__(self, *value):
        self.value = value
    def __str__(self):
        return repr(self.value)

