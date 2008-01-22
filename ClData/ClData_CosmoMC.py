"""
module/classes for CMB power spectrum data.

Initialize from COSMOMC data files -- see cosmomc/source/getdata.f90
"""

## to do

## *likelihood function for data given full C_l "theory" arrays
##   -or make this a model and define some sort of 'conversion' between the full C_l theory model
##    and other models such as Cosmo Params and binned C_l?

## *C_l theory class?

## should transpose Cl so it's [ClTT, ClTE, ClEE] (check order)

## *Some way to access the individual polarization spectra for a dataset
##   self.Clpol_idx[ch_type]=[start, end+1]  (need to enforce standard ch_type values...
##   -right now they're just in one contiguous array
##   -only *.newdat (B03) actual has polarization data now
##   -could retain this but just have accessor functions which put you in the right place.

from __future__ import division

import string
import os.path
import math
import operator
import copy

from numpy import (array, float64, int32, bool8, ones, zeros, nonzero, empty, rank,
                   reshape, log, exp, transpose, fabs, dot, arange, any, 
                   logical_not, where, logical_and)
import numpy.linalg as la

### "static" variables for this whole module
num_cls = 3   ### TT, TE, EE (change to 4 for BB)
lmax = 3100
halfsteps = 5

def initNumericalMarge(halfsteps):
    """ initialize numerical integration """
    step = arange(-halfsteps, halfsteps+1, 1.0, dtype=float64)
    margeweights = exp(-(step*3/float(halfsteps))**2/2)
    margenorm = sum(margeweights)
    return margenorm, margeweights

margenorm, margeweights = initNumericalMarge(halfsteps)  # make these 'static'        
class ClData_CosmoMC(object):
    """
    Data in CosmoMC format
    """
    
    ## what are the possibilities:
    ##   beam_uncertainty != 0        -- special-case for ==0?
    ##   calib_uncertainty != 0       --       ''
    ##   all_l_exact or has_xfactors  -- subclass based on this?
    ##       maybe not: the diffs are very localized?
    ##   has_corr_errors              -- special case for diag?
    ##                                   or subclass for quadratic forms?
    ##   has_pol

    ## what do we need to 'expose' to the world?
    ##   (Zth - Zobs) M^-1 (Zth - Zobs)
    ## where Z is either C_B or ln(C_B + x_B)
    ##    and Cth_B = sum_{l in B} Cth_l W_l
    ## so: *data* methods get_bandpowers(Cth)
    ##                    get_chi2(CB) (or lnLike?)
    ##      where the latter depends on has_xfactors
    ##        (or is subclassed depending on it?)
    ##        also need WMAP special case

    ## or maybe separate 'experiment' (which contains window fns)
    ##    from 'data'

    ## keep C_l, sig_l around even when xfactors are present for plotting purposes
    ## (look for use of 'copy' below)

    def __init__(self, file=None):
        if file is not None:
            self.readDataset(file)

    def readDataset(self, file):
        paramDict = readparams(file)
        self.name = paramDict['name']
        self.has_pol = str2bool(paramDict.get('has_pol', False))
        self.all_l_exact = str2bool(paramDict.get('all_l_exact', False))
        if self.all_l_exact:
            raise MyErr('All_l_exact not implemented')
            readAllExact()
            return
        self.num_points = int(paramDict.get('num_points', 0))
        self.calib_uncertainty = float(paramDict.get('calib_uncertainty', 0.0))
        self.beam_uncertain = str2bool(paramDict.get('beam_uncertainty', False))

        ## none of these needed outside of the setup
        window_dir = paramDict.get('window_dir')
           ### determine the absolute path?
           ### window_dir = "/".join([os.path.dirname(file), window_dir])
        windows_are_bare = str2bool(paramDict.get('windows_are_bare', False))
        windows_are_bandpowers = str2bool(paramDict.get('windows_are_bandpowers', True))
        windows_are_normalized = str2bool(paramDict.get('windows_are_normalized', False))

        file_points = int(paramDict.get('file_points', self.num_points))
        first_band = int(paramDict.get('first_band', 1))
        if first_band + self.num_points > file_points+1:
            raise ClDataError("first_band=", first_band,
                              " +num_points=", self.num_points,
                              " > file_points=", file_points)

        first_band -= 1   ### arrays are zero-index, so 0 is first possible

        if self.has_pol:
            self.ncls = num_cls
            self.Clpol_idx = {}
        else:
            self.ncls = 1
            self.Clpol_idx = {'TT': [0, self.num_points]}


            #print '%s: ncls: %d, first_band: %d, has_pol: %d, num_points: %d, beam_un: %d, calib_un: %d' % ( \
            #    self.name, self.ncls, first_band, \
            #    self.has_pol, self.num_points, self.beam_uncertain, self.calib_uncertainty)
            
        self.obs = empty(dtype=float64, shape=self.num_points)
        self.err_min = empty(dtype=float64, shape=self.num_points)
        self.err_pls = empty(dtype=float64, shape=self.num_points)
        self.beam_err = empty(dtype=float64, shape=self.num_points)
        self.window = zeros(dtype=float64,
                             shape=(self.num_points, self.ncls, lmax+1))
        self.win_min = empty(dtype=int32, shape=self.num_points)
        self.win_max = empty(dtype=int32, shape=self.num_points)
        self.ell = empty(dtype=int32, shape=self.num_points)

        ## or a full record at each point?
        
        band_file=str2bool(paramDict.get('bandpowers', False))
        if band_file:
            raise NotImplementedError
        else:
            for i in xrange(self.num_points):
                
                iline = i + first_band
                line = tuple([float(e) for e in
                              string.split(paramDict.get('data'+str(i+1)))])
                if self.beam_uncertain:
                    (self.obs[i], self.err_min[i], self.err_pls[i],
                     self.beam_err[i])= line
                else:
                    self.obs[i], self.err_min[i], self.err_pls[i] = line
                    self.beam_err[i]=0.0

                # print 'getting window array[%d], file[%d]' % (i, iline)
                self.readWindow(window_dir, i, iline, 
                                windows_are_bare,
                                windows_are_bandpowers,
                                windows_are_normalized)
                self.ell[i] = self.lmean(i)

        if self.beam_uncertain:
            self.beam_err /= self.obs

        self.sig = (self.err_min + self.err_pls)/2
        self.var = self.sig**2
        self.Cl = self.obs.copy()  #actual copy! even when has_xfactors
             ## save this in simple form in case xfactors are present
             ## could get fancy and use getattr to save memory?

        Ninv_file = paramDict.get('N_inv', None)
        self.has_corr_errors = Ninv_file is not None
        if self.has_corr_errors:
            tmp_mat = readMatrix(Ninv_file, file_points, file_points)
            if self.num_points != file_points:
                # truncate inverse
                tmp_mat = la.inv(tmp_mat)
                tmp_mat=tmp_mat[first_band:first_band+self.num_points,
                                first_band:first_band+self.num_points]
                self.N_inv = la.inverse(tmp_mat)
            else:
                self.N_inv = tmp_mat[0:self.num_points, 0:self.num_points]

        xfact_file = paramDict.get('xfactors', None)
        self.has_xfactors = xfact_file is not None
        if self.has_xfactors:
            ## data is ln(C+x), variance is sig^2/(C+x)^2
            tmp_x = readMatrix(xfact_file, self.num_points+first_band)
            self.has_xfactor = ones(dtype=bool8, shape=self.num_points)
            self.xfactors = tmp_x[first_band:first_band+self.num_points]
          #  print "Shapes -- xfactors:", self.xfactors.shape, \
          #        "  tmp_x: ", tmp_x.shape, "  var: ", self.var.shape, \
          #        "  obs: ", self.obs.shape
            self.var /= (self.obs+self.xfactors)**2
            self.obs = log(self.obs + self.xfactors)
            

    def readWindow(self, dir, i, iline, are_bare, are_bandpowers, are_normalized):
        """
        read the window functions in the directory dir
        nb. i = 0...num_points [index for array storage],
            iline = ifirst...ifirst+num_points [index in original bands]
        """

        win = self.window[i]
        fp = open("".join([dir,'/',self.name,str(iline+1)]))
        for lines in fp:
            splitline = lines.split()
            ll = splitline[0]
            l = int(float(ll))
            if fabs(l - float(ll)) > 1.0e-4:
                raise Nonint32Ell(i, ll)

            if l>=2 and l<=lmax:
                win[:,l] = [ float(w) for w in splitline[1:self.ncls+1] ]

        if l>lmax:
            print "   [ignored l=%d>%d]" % (l, lmax)

            
        if not are_bare:
            win[:,:] *= range(lmax+1)

        ## get the window min and max
        nz = where(win)[1]
        self.win_min[i], self.win_max[i] = min(nz), max(nz)

        if are_bandpowers:
            ###  include prefactor for sum_l C_l [W_l (2l+1)/4pi]
            
            ellwin = arange(self.win_min[i], self.win_max[i]+1)
            #for l in ellwin:
            #    win[:,l] *= (l+0.5)
            win[:, self.win_min[i]:self.win_max[i]+1] *= (ellwin+0.5)

            if not are_normalized:
                ## normalize to 1 = sum_l W_l (l+1/2)/(l(l+1)) 
                IW = sum(win[0,self.win_min[i]:self.win_max[i]+1]/(ellwin*(ellwin+1.0)))
                win[0,self.win_min[i]:self.win_max[i]+1] /= IW

            win[0,self.win_min[i]:self.win_max[i]+1] /= (2*math.pi) 

        if self.has_pol:
            self.inc_pol[i] = any(win[1,:] != 0)

        fp.close()
        
        ## remove small values
        smallfac = 1.0e-6
        maxw = max(fabs(win[0,:]))
        midx = where(logical_and(fabs(win[0,:])<maxw*smallfac, fabs(win[0,:]>0.0)))
        if len(midx[0])>0: 
            print "zeroing %d points< %f in window" % (len(midx[0]), maxw*smallfac)
            win[0,midx] *= 0

        # these methods all relate to "models" but require access to lots of local data

    def calcChisq(self, Cl):
        bandpowers = [self.getWinBandpower(i, Cl)
                      for i in range(self.num_points)]
        return self.getchisq(bandpowers)

    def calcLnLike(self,Cl):
        """
        compute -ln(likelihood(data | Cl));  (yes, negative)
        do numerical or analytic marginalization over beam and calibration uncertainty
        """

        if self.all_l_exact:
            raise MyErr('All_l_exact not implemented')
        else:
            denom = 1 
            bandpowers = array([self.getWinBandpower(i, Cl) for i in range(self.num_points)])

            if self.has_xfactors and (self.calib_uncertainty>1.0e-4 or self.beam_uncertain):
                chisq = self.getCalibMargexChisq(bandpowers)
            else:
                diffs = self.getdelta(bandpowers)
                chisq = self.getquadform(diffs)

                if self.calib_uncertainty > 1.0e-4 or self.beam_uncertain:
                    tmp = self.getlinform(bandpowers)
                    chi2op = dot(diffs, tmp)
                    chi2pp = dot(bandpowers, tmp)
                    if self.beam_uncertain:
                        beam = self.beam_err*bandpowers
                        tmp = self.getlinform(beam)
                        chi2dd = dot(beam, tmp)
                        chi2pd = dot(bandpowers, tmp)
                        chi2od = dot(diffs, tmp)

                    if self.calib_uncertainty > 1.0e-4:
                        # analytic marginalization over calib. uncertainty
                        wpp = 1/(chi2pp+1/self.calib_uncertainty**2)
                        chisq -= wpp*chi2op**2
                        denom /= wpp*self.calib_uncertainty**2
                    else:
                        wpp = 0
                        
                    if self.beam_uncertain:
                        wdd = 1/(chi2dd - wpp*chi2pd**2 + 1)
                        chisq -= wdd*(chi2od-wpp*chi2op*chi2pd)**2
                        denom /= wdd
                        
            if denom!=1:
                chisq += log(denom)

        return chisq/2

    def getCalibMargexChisq(self, bandpowers):
        """
        Numerically integrate over the calibration uncertainty
        Assume Gaussian prior, as for analytic calculation without x-factors
        """

        bandpowers = array(bandpowers)

        chisq=empty(dtype=float64, shape=2*halfsteps+1)
        chisqcalib=empty(dtype=float64, shape=2*halfsteps+1)
        low = -30.0 + zeros(2*halfsteps+1, float64)

        hrange=range(-halfsteps, halfsteps+1)
        for ibeam in hrange:

            beambandpowers = bandpowers
            if self.beam_uncertain:
                beambandpowers *= 1+self.beam_err*ibeam*3/halfsteps

            for i in hrange:
                calib = 1+self.calib_uncertainty*i*3/halfsteps #Go out to 3 sigma
                chisq[i+halfsteps] = self.getchisq(calib*beambandpowers)

            minchisq = min(chisq)
            #### deal with underflow
            exparg = array([max(z) for z in zip(low, -(chisq-minchisq)/2)])
            chisqcalib[ibeam+halfsteps] = -2*log(sum( \
                margeweights*exp(exparg)/margenorm)) + minchisq

            if not self.beam_uncertain:
                return chisqcalib[ibeam+halfsteps]

        minchisq = min(chisqcalib)
        exparg = array([max(z) for z in zip(low, -(chisqcalib-minchisq)/2)])
        return  -2*log(sum(margeweights*exp(exparg)/margenorm)) + minchisq


    ## modified to deal with BP+x<0
    def getdelta(self, BP):
        """
        get the appropriate (theory - data) difference
        depending on presence of xfactors
        """
        diffs=empty(dtype=float64, shape=self.num_points)
        if self.has_xfactors:
            idx = self.has_xfactor
            ndx = logical_not(idx)
            
            zth = BP[idx] + self.xfactors[idx]
            #zth[where(zth<=0.0)] = 1.0e-10   ### AHJ

            diffs[idx] = self.obs[idx] - log(zth)
            diffs[ndx] = self.obs[ndx] - BP[ndx]
        else:
            diffs = self.obs - BP

        return diffs

    def getchisq(self, BP):
        """
        get the chi^2 associated with this set of BP=bandpowers,
        adjusting for xfactors and correlated noise
        """
        return self.getquadform(self.getdelta(BP))

    def getlinform(self, vec):
        """ get the linear form vec^T.N_inv or sum(vec/sig^2) """
        if self.has_corr_errors:
            return dot(self.N_inv,vec)
        else:
            return vec/self.var
        

    def getquadform(self, vec):
        """ get the quadratic form vec^T.N_inv.vec or sum(vec^2/sig^2)"""
        if self.has_corr_errors:
            return dot(vec,dot(self.N_inv,vec))
        else:
            return sum(vec**2/self.var)

    def getWinBandpower(self, i, Cl):
        """
        return the bandpower if observing Cl for this experiment in band i
        nb Cl = [ClTT, ClTE, ClEE], window[i] = [WTT, WTE, WEE]
        and i just numbers band; doesn't distinguish between TT/TE/EE
        """

        if rank(Cl)==2:
            ClTT = Cl[0]
        elif rank(Cl)==1:
            ClTT=Cl
        else:
            raise MyErr('Cl shape wrong in getWinBandpower:', Cl.shape)

        maxl = min([len(ClTT), self.win_max[i]+1])
        bandpower = dot(ClTT[self.win_min[i]:maxl],
                        self.window[i, 0, self.win_min[i]:maxl])
        if self.has_pol and rank(Cl)==2 and self.inc_pol[i]:
            bandpower += (Cl[1:num_cls, win_min[i], maxl+1]*
                          self.window[i, 1:num_cls, win_min[i], maxl+1]).sum()
            print 'Accessing polarization getWinBandpower: band', i
            #for l in xrange(self.win_min[i], maxl):
            #    bandpower += dot(Cl[1:num_cls, l], self.window[i, 1:num_cls, l])            
        return bandpower

    def lmean(self, band):
        """ average <ell> for the band
            need to correct for polarized channels
            -- only include the appropriate window
        """
        return dot(range(self.win_min[band], self.win_max[band]+1),
                   self.window[band, 0, self.win_min[band]:
                               self.win_max[band]+1])/ \
                   sum(self.window[band, 0, self.win_min[band]:
                               self.win_max[band]+1])


### module-level functions from here down..............

def readparams(file=None):
    """
    read a parameter with lines of the form
       param = a very long value with spaces allowed
    return a dict with dict[param]=value
    empty lines and lines beginning with ; or # are ignored
    convert key name to lowercase
    """
    fp = open(file)  # or just allow already opened file-like-object?

    dc = {}
    for line in fp:
        key, val = parseline(line)
        if key is not None:
            dc[string.lower(key)] = val
    fp.close()
    return dc


def parseline(line):
    """
    parse a line of the form
       param = Value With Whitespace
    into the tuple ('param', 'Value With Whitespace')
    return None, None if blank line or first non-blanck char is ';' or '#'
    raise ParamSyntaxError if no '='
    """
    line = line.strip()
    if len(line) == 0 or line[0] == ';' or line[0] == '#':
        return None, None
    splt = string.split(line, '=', 1)
    if len(splt) != 2:
        raise ParamSyntaxError(line)
    return splt[0].strip(), splt[1].strip()

def str2bool(str):
    """ convert a string like 'T', 'True', 'true' etc to boolean
        just return if already boolean
        very little error checking!"""
    if str in [True, False]: return str
    return string.upper(str[0])=='T'

def readMatrix(filename, *shape):
    """
    read a matrix of given shape [e.g., (n,m) for an n*m matrix]
        stored as ascii in filename
    special case: leave as-is if no shape given
    (should probably work on file pointers too?)
    if one-dimensional, just keep the first n values
    """
    f = open (filename)
    data=[]
    for line in f:
        linedata = []
        for val in line.split():
            linedata.append(float(val))
        if len(linedata):
            data.append(linedata)
    f.close()
    if shape is ():
        return array(data)
    else:
        if len(shape)==1:
            return reshape(array(data)[0:shape[0]], shape)
        elif len(data) == int(reduce(operator.mul, shape, 1.0)):
            return reshape(array(data), shape)
        else:
            raise ArraySizeError(filename, shape, data)
        
            

class MyErr(Exception):
    def __init__(self, *value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class ParamSyntaxError(MyErr):
    pass

class ClDataError(MyErr):
    pass

class Nonint32Ell(MyErr):
    pass

class ArraySizeError(MyErr):
    pass

