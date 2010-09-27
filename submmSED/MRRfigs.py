#!/usr/bin/env python
# encoding: utf-8
"""
MRRfigs.py

Created by Andrew H. Jaffe on 2010-09-25.
Copyright (c) 2010 Imperial College London. All rights reserved.
"""

import sys
import os

import cPickle as pickle
import driver

import matplotlib.pyplot as plt
import numpy as np

def figs(ret123=None, mean_or_ML='mean',lab=""):
    
    # with open("out_MRR/out_[0].pickle") as f: r0=pickle.load(f)[0]
    # with open("out_MRR/out_[1].pickle") as f: r1=pickle.load(f)[1]
    # with open("out_MRR/out_[2].pickle") as f: r2=pickle.load(f)[2]
    
    if ret123 is None:
        ret123 = driver.postprocess("./out_MRR/")
    
    ifig = 0
        
    suff = lab+".png"
        
    plt.figure(ifig)
    T = ret123[1][:][mean_or_ML][:][:,2]
    Tmean = np.mean(T)
    Tsig = np.std(T)
    lab = r"$T = %5.2f \pm %5.2f$" % (Tmean, Tsig)
    plt.hist(T, bins=20)
    plt.xlabel("Temperature (K)")
    plt.title(lab)
    plt.savefig("THist"+suff)

    ifig +=1; plt.figure(ifig)
    beta = ret123[1][:][mean_or_ML][:][:,1]
    bmean = np.mean(beta)
    bsig = np.std(beta)
    lab = r"$\beta = %5.2f \pm %5.2f$" % (bmean, bsig)
    plt.hist(beta, bins=20)
    plt.xlabel(r"$\beta$")
    plt.title(lab)
    plt.savefig("betaHist"+suff)
    
    ifig+=1; plt.figure(ifig)
    Tb2 = ret123[2][:][mean_or_ML][:][:,-1]
    Tb2mean = np.mean(Tb2)
    Tb2sig = np.std(Tb2)
    lab =  r"$T = %5.2f \pm %5.2f$ (fixed $\beta=2$)" % (Tb2mean, Tb2sig)
    plt.hist(Tb2, bins=20)
    plt.xlabel("Temperature (K)")
    plt.title(lab)
    plt.savefig("THist_beta2"+suff)
        
    ifig+=1; plt.figure(ifig)
    plt.plot(T, Tb2, ',')
    plt.xlabel("Temperature (K)")
    plt.ylabel(r"Temperature (K) [fixed $\beta=2$]")
    plt.plot([0,45], [0,45])
    plt.title(lab)
    plt.savefig("SingleT_betafix"+suff)
    
    ifig+=1; plt.figure(ifig)
    plt.plot(T, beta, ',')
    plt.xlabel("Temperature (K)")
    plt.ylabel(r"$\beta$")
    plt.savefig("Tbeta"+suff)

    ifig+=1; plt.figure(ifig)
#    objs = np.logical_and(ret123[0][:][mean_or_ML][:][:,0] > -10, ret123[0][:][mean_or_ML][:][:,2] > -10)
    ### below is not very sensitive to the limit used
    objs = (ret123[0][:][mean_or_ML][:][:,0] - ret123[0][:][mean_or_ML][:][:,2] > -20)
    print "got %d 2 T objs" % sum(objs)
    T1 = ret123[0][objs][mean_or_ML][:][:,1]
    T2 = ret123[0][objs][mean_or_ML][:][:,3]
    T1mean = np.mean(T1)
    T1sig = np.std(T1)
    T2mean = np.mean(T2)
    T2sig = np.std(T2)
    lab = "$T_1 = %5.2f \pm %5.2f$, $T_2 = %5.2f \pm %5.2f$" % (T1mean, T1sig, T2mean, T2sig)
    
    plt.plot(T1, T2, ',')
    plt.xlabel("Temperature (K)")
    plt.ylabel(r"Temperature (K)")
    plt.title(lab)
    plt.plot([0,30], [0,30])
    plt.savefig("TwoT"+suff)
    
    lowTobjs = np.logical_and(objs, ret123[0][:][mean_or_ML][:][:,1]<10)
    lowTnames = ret123[0][lowTobjs]['name'][:]
    with open("lowT.txt", 'w') as f:
        for n in lowTnames:
            f.write(n + "\n")
            
            
    return ret123
    
    