#!/usr/bin/env python
# encoding: utf-8
"""
Clfigs.py

Created by Andrew Jaffe on 2008-05-03.
Copyright (c) 2008 __MyCompanyName__. All rights reserved.
"""

from __future__ import with_statement

import sys
import os
import numpy as N
import pylab
import matplotlib.pyplot as plt


def figs(rel=True, single=False, wmap=False, ylim=10000):
    
    wmapfile = "/Users/jaffe/Desktop/Projects/GUS/wmap_binned_tt_spectrum_5yr_v3.txt"
    
    wm = N.loadtxt(wmapfile, unpack=True)
    
    gdir = "../dat"
    prefix = ["bolo/bolo", "interf/interf", "radiometer/rad", "wang1/wang2","ground.2/ground"]
    labels = ['Bolometers','Interferometers','Radiometers','Wang et al', 'sub-orbital']
    
    tmp = N.fromfile(gdir+"/" + "CarloClModel.dat", sep=" ")
    tmp.shape = -1, 6
    ells = N.int_(tmp.T[0])  ## nb doesn't usually start at l=0
    
    max_ell = N.max(ells)
    ll = N.arange(max_ell+1)
    llClTT = N.zeros(max_ell+1)
    llClEE = N.zeros(max_ell+1)
    llClTE = N.zeros(max_ell+1)
    llClTT[ells] = tmp.T[1]
    llClEE[ells] = tmp.T[2]
    llClTE[ells] = tmp.T[3]
    
    
    sym =    ".xo+s^"
    colors = "bgrcmy"
    l = []; C = []; e = []; chi2 = []
    for i, pre in enumerate(prefix):
        if single:
            plt.figure()
            if not rel:
                plt.plot(ll, llClTT)
        with open(gdir+"/"+pre+".bp") as f:
            dat = N.loadtxt(f,unpack=True)
            l.append(N.array(dat[0], dtype=int))
            C.append(dat[1])
            e.append(dat[2])
            if rel:
                C[-1][:] /= llClTT[l[-1]]
                e[-1][:] /= llClTT[l[-1]]
                chi2.append((((C[-1]-1)/e[-1])**2).mean())
                print labels[i]+":", chi2[-1]
                lab = "%s [%4.2f]" % (labels[i],chi2[-1])
            else:
                lab = labels[i]
            plt.plot(l[-1],C[-1],colors[i]+sym[i], label=lab)
            plt.errorbar(l[-1],C[-1],e[-1],fmt=colors[i]+sym[i])
            if single:
                plt.plot(wm[0], wm[3], 'xk', label='WMAP')
                plt.errorbar(wm[0], wm[3], wm[4], fmt='xk')
            plt.xlim(0,2500)
            if ylim: plt.ylim(0,10000)
            if single:
                plt.legend()
    if rel:
        plt.axhline(1)
        plt.ylim(0,4)
        plt.legend()
    elif not single:
        if wmap:
            plt.plot(wm[0], wm[3], 'xk', label='WMAP')
            plt.errorbar(wm[0], wm[3], wm[4], fmt='xk')
        plt.plot(ll, llClTT)
        plt.legend()
        plt.xlim(0,2500)
        if ylim: plt.ylim(0,10000)


def main():
	pass


if __name__ == '__main__':
	main()

