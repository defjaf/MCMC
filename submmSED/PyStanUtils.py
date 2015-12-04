"""

utilities for working with pystan and interfacing with my MCMC code and AL's getdist

should split this file into generic MCMC/pystan/getdist functions and very specific submm ones

"""

from __future__ import print_function

import warnings
import csv
import sys
import os
import itertools

## python 2/3 compatibility
try:
    import cPickle as pickle
except ImportError:
    import pickle
import gzip
try:
    from importlib import reload
except ImportError:
    pass

from six import iteritems

import matplotlib.pyplot as plt
import numpy as np
import pystan
import getdist
import getdist.plots

sys.path.append(os.environ['HOME']+'/home/proj/stats/MCMC')

import model     ## needed for plotting SED models (duplicates some STAN code)

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
    

## Plot pystan output using MCMC/subMM routines
def plot_pystan(data, fit, ncomp=None, model_name=None, linear=True, wavelength=True, logplot=True, label=None):
    """ 
    plot the output from a pystan run as a spectrum
       nb. can combine several calls into one figure.
    data: an actual AHJ MCMC data object (needs to have .plot method)
    fit: pystan fit object
    model_name: name of pystan model (from fit if None)
    ncomp: number of components (1 or 2; should be 1 for optically thick model)
           from fit if None
    label: label for plot legend. Need to call plt.legend() explicitly
    """
    
    params = fit.get_posterior_mean().mean(axis=1)[:-1]
    flatnames = fit.flatnames
    
    if model_name is None:
        model_name = fit.model_name
        model_name = model_name[0:model_name.rfind("_")]
        
    if ncomp is None:
        ncomp = fit.data["N_comp"]

    if model_name in ("greybody", "greybody_beta2"):
        assert ncomp in (1,2)
        paramorder = ["amplitude[0]", "beta[0]", "T[0]"]
        if ncomp == 1:
            mod = model.submmModel1_normalized_logA if not linear else model.submmModel1_normalized
        elif ncomp == 2:
            mod = model.submmModel2_normalized_logA if not linear else model.submmModel2_normalized
            paramorder += ["amplitude[1]", "beta[1]", "T[1]"]
    elif model_name == "greybody_thick":
        assert ncomp==1
        mod = model.submmModel1_opticallythick_logA  if not linear else model.submmModel1_opticallythick
        paramorder = ["amplitude[0]", "beta[0]", "T[0]", "nu0[0]"]

    if model_name == "greybody_beta2":
        ### add in beta=2 fixed parameters
        params = params.tolist()
        for i in range(ncomp):
            flatnames.append("beta[%d]" % i)
            params.append(2.0)
        params = np.array(params)

    assert len(params)==len(paramorder)
    assert len(params)==len(flatnames)
    flatnames = np.array(flatnames)

    ### package the params for my MCMC model code (duplicates stuff done in stan code...)
    MPparams = [params[flatnames==p] for p in paramorder]

    MPmod = mod(*MPparams)
    MPmod.plot(data, wavelength=wavelength, logplot=logplot, label=label)
        
        
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

def pystan_postprocess_togetdist(allfits, texdict=texdict, do_pickle=False):
    """ convert pytstan output to getdist  
        -- generic, requires output in dictionary """

    gds = {}

    for obj, fs in iteritems(allfits):
        gds[obj] = [pystan2getdist(f, texdict) for f in fs]
    ### order is 1comp, 2comp, 1compb2, 2compb2, thick

    if do_pickle:
        with gzip.open("greybody_getdist.pkl.gz", "wb") as fil:
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
                
                
def pystan_postprocess_SED_plot(allfits, label=""):
    """ plot pytstan submm sed output as SEDs  -- specific to greybody models """

    labs = ("1 comp", "2 comp", r"1 $\beta=2$", r"2 $\beta=2$", "thick")
    for objname, fits in iteritems(allfits):
        for fit, lab in zip(fits, labs):
            plot_pystan(dat[objname], fit, label=lab)
            plt.title(objname)
        plt.legend(loc='best')

        plt.savefig("%s_%sSED.png" % (objname, label), dpi=200)
        plt.figure()
        
        
        
def pystan_postprocess_csv(allfits, dat, label="DLC_"):
    """ convert pystan submm sed output to csv -- specific to greybody models """

    obj_names = dat.keys()

    modnames = ["1comp","2comp","1compb2","2compb2","thick"]
    for i, mod in enumerate(modnames):
        f = [a[i] for a in allfits.values()]
        fname = label + mod + ".csv"
        ptab = pystan2table(f)
        table2csv(fname, ptab, obj_names)