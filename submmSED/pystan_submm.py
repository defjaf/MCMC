
from __future__ import print_function

import sys
import os

## python 2/3 compatibility
from six import iteritems

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.environ['HOME']+'/home/proj/stats/MCMC')

import pystan_utils

import model     ## needed for plotting SED models (duplicates some STAN code)

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
    
    
def pystan_postprocess_SED_plot(allfits, dat, label=""):
    """ plot pytstan submm sed output as SEDs  -- specific to greybody models 
    assumes each entry in allfits is in the order of `labs` for the legend
    TODO: probably need to be able to specify a different model order
    """

    labs = ("1 comp", "2 comp", r"1 $\beta=2$", r"2 $\beta=2$", "thick")
    for objname, fits in iteritems(allfits):
        for fit, lab in zip(fits, labs):
            plot_pystan(dat[objname], fit, label=lab)
            plt.title(objname)
        plt.legend(loc='best')

        plt.savefig("%s_%sSED.png" % (objname, label), dpi=200)
        plt.figure()
        
        
        
def pystan_postprocess_csv(allfits, dat, label="DLC_", modnames=None):
    """ convert pystan submm sed output to csv -- specific to greybody models 
    currently assumes that there are entries for each of the models. 
    
    TODO: allow not having all possible models. 
    see *_SED_plot() above -- 
    as in that case, we probably need to be able to specify a different model order
    """

    obj_names = dat.keys()

    if modnames is None:
        modnames = ["1comp","2comp","1compb2","2compb2","thick"]
    nmod = len(allfits.values()[1])  ### assumes same number of models for each object
    for i in range(nmod):
        f = [a[i] for a in allfits.values()]
        if nmod == len(modnames): 
            mod = modnames[i]
        else:
            mod = str(i)
        fname = label + '_' + mod + ".csv"
        ptab = pystan_utils.pystan2table(f)
        pystan_utils.table2csv(fname, ptab, obj_names)