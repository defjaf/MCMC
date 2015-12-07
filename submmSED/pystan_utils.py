"""

utilities for working with pystan and interfacing with my MCMC code and AL's getdist

should split this file into generic MCMC/pystan/getdist functions and very specific submm ones

"""

from __future__ import print_function

import csv
import itertools

## python 2/3 compatibility
try:
    import cPickle as pickle
except ImportError:
    import pickle
import gzip

from six import iteritems

import numpy as np
import getdist
import getdist.plots

# convert from my MCMC/subMM data format to pystan
def makeStanDict(dobj, N_comp=2):
    """ convert from an MCMC data object to a dict for input to STAN"""
    return dict(
        N_comp = N_comp,
        N_band = dobj.n,
        nu_obs = dobj.freq_obs,
        flux = dobj.d,
        sigma = dobj.sig,
        z = dobj.z
    )
    

### pystan output to A Lewis getdist
texdict={"amplitude": "A", "beta": "\\beta", "tau": "\\tau", "nu0": "\\nu_0"}
def pystan2getdist(fit, texdict=texdict):
    """ convert a pystan fit object into an A Lewis MCSamples object
        texdict is a dictionary mapping from stan parameter names to latex labels, if desired
    """
    
    ### have to either deal with combining chains (permute=False) 
    ###                       or separating vector-valued parameters (permute=True)
    
    sample_dict = fit.extract(permuted=True)
    names = []
    labels = []
    samples = []
    for par, samps in iteritems(sample_dict):
        if par == "lp__":
            loglikes = -samps
        else:
            if len(samps.shape)==1:
                names.append(par)
                labels.append(texdict.get(par,par))
                samples.append(samps)
            else:
                for i,s in enumerate(samps.T):
                    samples.append(s)
                    names.append(par+str(i+1))
                    
                    labels.append(texdict.get(par,par)+"_"+str(i+1))
    
    return getdist.MCSamples(names=labels, samples=samples, loglikes=loglikes, labels=labels)
    

        
        
##  convert pystan text and param output

def pystan2table(fits):
    """ convert list of pystan fit object into a parameter table -- generic """

    parnames = fits[0].flatnames + ["lp__",]

    npar = len(parnames)
    nobj = len(fits)

    param_table = np.empty((nobj,), dtype={"names": parnames, "formats": [np.float,]*npar})
    
    for i, fit1 in enumerate(fits):
        param_table[i] = fit1.get_posterior_mean().mean(axis=1)
        
    return param_table
    
    
    
##names = dat.keys()
def table2csv(fil, table, names=None):
    """ convert a structured numpy array into a csv, possibly prepending a 'names' column """
    names = [None,] if not names else names
    with open(fil, 'wb') as f:
        fw = csv.writer(f)
        fw.writerow(("object",) + table.dtype.names )
        for obj, row in itertools.izip_longest(names, table):
            fw.writerow( [obj,] + [r for r in row] )
            
            
### postprocessing
    
def pystan_postprocess_text(allfits, outfile):
    """ text output from pystan 
        -- generic, requires output in dictionary """ 
    with open(outfile, "w") as f:
        for name, fit_1obj in iteritems(allfits):
            print("============="+name+"=============", file=f)
            for i, fit_1model in enumerate(fit_1obj):
                try:
                    print(fit_1model, file=f)
                except OverflowError:
                    print("***** Can't print: %d %s *****" % (i, name), file=f)
                    print(fit_1model.summary())
                    
                    
                    
texdict={"amplitude": "A", "beta": "\\beta", "tau": "\\tau", "nu0": "\\nu_0"}

def pystan_postprocess_togetdist(allfits, texdict=texdict, 
                                 do_pickle=False, picklename="greybody_getdist.pkl.gz"):
    """ convert pytstan output to getdist  
        -- generic, requires output in dictionary with a sequence at each entry """

    gds = {}

    for obj, fs in iteritems(allfits):
        gds[obj] = [pystan2getdist(f, texdict) for f in fs]
    ### order is 1comp, 2comp, 1compb2, 2compb2, thick

    if do_pickle:
        with gzip.open(picklename, "wb") as fil:
            pickle.dump(gds, fil)
    return gds
    
    
def pystan_postprocess_getdist_plot(gds, label=""):
    """ getdist plots from pystan (after conversion)  
        -- generic, requires output in dictionary """

    for obj, gd in iteritems(gds):
        print(obj)
        for i, gd1 in enumerate(gd):
            g = getdist.plots.getSubplotPlotter()
            try:
                g.triangle_plot(gd1)
                g.export("%s_%d%s.png" % (obj, i, label))

            except Exception as E:
                print("at %s %d, get exception" % (obj, i))
                print(E)
                
                