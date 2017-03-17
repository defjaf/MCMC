from __future__ import division

import os.path
import math
import string
import copy

from numpy import (array, exp, log, transpose, zeros, ones, 
                   int32, float64, bool8, empty, concatenate)
import numpy.linalg as la

import ClData_CosmoMC
from ClData_CosmoMC import num_cls, lmax

class ClData_bcp(ClData_CosmoMC.ClData_CosmoMC):
    """
    represents *.newdat CMB data a la bcp
    """
    
    def readDataset(self, file):
        """
        Read a dataset in the following format from Carlo Contaldi
        (courtesy CosmoMC)
        !name
        !n_bands_TT n_EE, n_BB, n_EB, n_TE, n_TB
        !has_calib_uncertain calib calib_err
        !has_beam_uncertain beam beam_err
        !ilike (0: Gaussian, 1: all x-factor, 2: specified have x-factor)
        !loop over {
        ! band-types
        ! band info: num obs + - x_factor l_min l_max use_x
        ! correlation-matrix (ignored)
        ! }  
        ! covariance matrix
        """

        minmax = empty(dtype=int32, shape=(2,6))
        
        fp = open(file)  # need dir information?
        self.name = fp.readline().strip()
        fisher_T_cmb = False
        if self.name == 'FISHER_T_CMB':
            fisher_T_cmb=True
            self.name = fp.readline().strip()

        line = string.split(fp.readline())  # number of points in TT EE BB ?? TE ??
        npol = [int(nstr) for nstr in line]
        self.has_pol = (sum(npol[1:]) > 0)

        self.all_l_exact=False
        windows_are_bare = False
        windows_are_bandpowers = True
        windows_are_normalized = True
        
        file_points = sum(npol)
        self.num_points = sum(npol)
        n_types = sum([1 for nn in npol if nn>0])  ## number of nonzeros
        line = fp.readline().strip()
        if line == 'BAND_SELECTION':
            #list of 'first_band last_band' for each pol type
            #if first_band=0 then ignore that pol type
            self.num_points = 0
            self.has_pol = False
            for i in range(6):
                minmax[:,i] = [int(mm) for mm in fp.readline().split()]
                if minmax[0,i] != 0:
                    self.num_points += minmax[1,i] - minmax[0,i] + 1
                    self.has_pol = (i > 0)
                else:
                    minmax[1,i] = 0
            
            line = fp.readline().strip() 
        else:
            #use all bands in file
            for i in range(6):
                minmax[:,i] = [1, npol[i]]

        minmax -= 1  ## convert to 0-based indexing

        cal, self.calib_uncertainty = [ float(e) for e in line.split()[1:3] ]

        line = fp.readline().split()
        beam_width, beam_sigma =  [ float(e) for e in line[1:3] ]
        self.beam_uncertain = (int(line[0]) != 0)

        if self.has_pol:
            self.ncls = num_cls
        else:
            self.ncls = 1
        self.obs = empty(dtype=float64, shape=self.num_points)
        self.err_min = empty(dtype=float64, shape=self.num_points)
        self.err_pls = empty(dtype=float64, shape=self.num_points)
        self.beam_err = empty(dtype=float64, shape=self.num_points)
        self.window = zeros(dtype=float64,
                             shape=(self.num_points, self.ncls, lmax+1))
        self.inc_pol = empty(dtype=bool8, shape=self.num_points)
        self.win_min = empty(dtype=int32, shape=self.num_points)
        self.win_max = empty(dtype=int32, shape=self.num_points)
        self.ell = empty(dtype=int32, shape=self.num_points)
        used_bands = []
        tmp_x = empty(dtype=float64, shape=self.num_points)
        lb = empty(dtype=float64, shape=(2, self.num_points))
                   
        ilike = int(fp.readline().split()[0])
        self.has_xfactors = ilike>0
        #  print 'ilike=%d, has_xfactors=%d' % (ilike, self.has_xfactors)
        #  1 : all bands are offset lognormal
        #  2 : only bands specified have offset lognormal
        if ilike>0:
            self.has_xfactor = ones(dtype=bool8, shape=self.num_points)

        use_i = 0   ## keeps track of the global band index = range(self.num_points)
        file_i = 0

        self.Clpol_idx = {}
        
        inpol = [k for k in xrange(6) if npol[k]!=0]
        ## maps from all 6 polz types to the ones actually present
        for k in inpol:
            ch_type = fp.readline().strip()[0:2]
            self.Clpol_idx[ch_type] = [use_i, use_i]
            print '  ch_type=%s' % ch_type
            for i in xrange(npol[k]):
                line = fp.readline().split()
                if i>=minmax[0,k] and i<=minmax[1,k]:
                    self.Clpol_idx[ch_type][1] += 1
                    used_bands.append(file_i)
                    (self.obs[use_i], self.err_min[use_i], self.err_pls[use_i], \
                     tmp_x[use_i], lb[0,use_i], lb[1,use_i] \
                     ) = [float(e) for e in line[1:7]]
                    if ilike>1:
                        self.has_xfactor[use_i] = int(line[7])

                    # print 'getting window array[%d], file[%d], global[%d])' % (use_i, file_i, i)
                    ## nb. the window functions are always W[iCl, l] for iCl=[TT,TE,EE[,BB]]
                    self.readWindow('data/windows', use_i, file_i, \
                                                    windows_are_bare, \
                                                    windows_are_bandpowers, \
                                                    windows_are_normalized)
                    self.ell[use_i] = sum(lb[:, use_i])/2.0
                    use_i += 1  ## careful of alignment here (inside 'if i')
                file_i += 1     ## (inside 'for i') 
                
            ## discard the correlation matrix
            for i in xrange(npol[k]):
                fp.readline()

        self.has_corr_errors = True
        self.N_inv = empty(dtype=float64, shape=(self.num_points, self.num_points))

        ### get the rest of the file into tmp_mat
        nignore = 0
        lines = fp.readlines()
        tmp_lis=[]
        for line in lines:
            try:
                tmp_lis.append([float(e) for e in line.split()])
            except ValueError:
                nignore += 1
                #print "Possible syntax error; ignoring line:", line.strip()
            
        if nignore:
            print "Ignored %d lines." % nignore

        tmp_lis = concatenate(tmp_lis)
        tmp_mat=array(tmp_lis).copy()
        tmp_mat.shape = (file_points, file_points)
        
        ##print tmp_mat.shape

        # is this efficient?
        idx = [ used_bands ] * self.num_points   ## replicate, not multiply!
        self.N_inv=tmp_mat[ transpose(idx), idx ]
        
   #     print "Noise mat, first, last entries:", \
   #        self.N_inv[0,0], self.N_inv[self.num_points-1, self.num_points-1]


        self.beam_err = exp(-self.ell*(self.ell+1)*
                            1.526e-8*2.0*beam_sigma*beam_width)-1.0
        self.beam_err = abs(self.beam_err)
        self.sig = (self.err_pls + self.err_min)/2.0
        self.obs *= cal**2
        self.Cl = self.obs.copy()  #actual copy! even when has_xfactors
        self.sig *= cal**2
        self.var = self.sig**2

        self.N_inv *= cal**4
        if fisher_T_cmb:
            self.N_inv *= 2.725**4 * 1.0e24

        if self.has_xfactors:
            self.xfactors = cal**2 * tmp_x[0:self.num_points]
            ## transform to z=ln(C+x) for N_inv,  obs and var

            ## there must be a more pythonic way to do this than these loops!
            ## (but obvious ways with index arrays don't work...
            for i in range(self.num_points):
                for j in range(self.num_points):
                    if self.has_xfactor[i]:
                        self.N_inv[i,j] /= self.obs[i]+self.xfactors[i]
                    if self.has_xfactor[j]:
                        self.N_inv[i,j] /= self.obs[j]+self.xfactors[j]
                            
            
            idx = self.has_xfactor            
            self.var[idx] /= (self.obs[idx]+self.xfactors[idx])**2
            self.obs[idx] = log(self.obs[idx] + self.xfactors[idx])
                    
        self.N_inv = la.inv(self.N_inv)
        
