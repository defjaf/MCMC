"""
read the tables produced by driver.writeTab

may require hacking the text files to make sure there is a space before the (z=*) entry

"""

import numpy as np
import matplotlib.pyplot as plt

def readtab(f):
    with open(f) as fp:
        hdr = fp.readline()
        
        ## header is names separated by >1 space
        names = [s.strip() for s in hdr.split("  ") if s]
        
        ret = np.genfromtxt(fp, names=names)
        
    for col in ret.dtype.fields.keys():
        print("%s = %f +- %f" % (col, ret[col].mean(), ret[col].std()))

    return ret
        
def plotpeel(f="./out_peel/0.npy", col="Mean_param_1", dat="./pixelfit.dat"):
    tab = readtab(f)
    i,j = np.loadtxt(dat, usecols=[0,1], unpack=True)
    plt.scatter(i,j,40,c=tab[col],marker='s')
    return i,j, tab

def getdat():

    n0 = np.genfromtxt("./submmSED/out_next0/next0_00.npy", skiprows=1, dtype=None, unpack=False)
    n1 = np.genfromtxt("./submmSED/out_next0/next0_01.npy", skiprows=1, dtype=None, unpack=False)

    d0 = np.genfromtxt("./submmSED/out_XX/fullcat_0.npy", skiprows=1, dtype=None, unpack=False)
    d1 = np.genfromtxt("./submmSED/out_XX/fullcat_1.npy", skiprows=1, dtype=None, unpack=False)

    n0names = [o[0] for o in n0]
    n1names = [o[0] for o in n1]

    d0 = [o for o in d0 if o[0] in n0names] ### assumes same order
    d1 = [o for o in d1 if o[0] in n1names] ### assumes same order

    n0 = np.array([list(o)[2:] for o in n0])
    n1 = np.array([list(o)[2:] for o in n1])
    d0 = np.array([list(o)[2:] for o in d0])
    d1 = np.array([list(o)[2:] for o in d1])
    
    return n0names, n1names, n0, n1, d0, d1

## model 0: columns 1-4 for ML, 5-8 for mean

## model 1: 1-3 for ML, 4-6 for mean

def plotcomp(n, d, cols):
    
    for i, c in enumerate(cols):
        
        plt.subplot(2,2,i+1)
        plt.plot(n[:,c], d[:,c], ',')
    
        
def main():
    n0names, n1names, n0, n1, d0, d1 = getdat()
    
    plt.figure(0)
    plotcomp(n0, d0, [1,2,3,4])
    
    plt.figure(1)
    plotcomp(n0, d0, [5,6,7,8])
    
    plt.figure(2)
    plotcomp(n1, d1, [1,2,3])
    
    plt.figure(3)
    plotcomp(n1, d1, [4,5,6])

