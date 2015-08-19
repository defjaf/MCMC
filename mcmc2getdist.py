"""
Convert an AHJ MCMC object into an A Lewis MCSamples object
"""

import getdist

def mcmc2getdist(mcmc):
    """
    Convert an AHJ MCMC object into an A Lewis MCSamples object
    """

    ## only use parameters with proposal width>0 -- otherwise fixed
    params = mcmc.prop.sigmas>0

    samples = mcmc.samples[:,params]
    labels = [t for t,p in zip(mcmc.mod.texNames, params) if p]
    loglikes = -mcmc.lnPr   ### note negative: AL version wants -lnPr
    return getdist.MCSamples(names=labels, samples=samples, loglikes=loglikes, labels=labels)

def pystan2getdist(fit):
    """ convert a pystan fit object into an A Lewis MCSamples object """
    
    ### have to either deal with combining chains (permute=False) 
    ###                       or separating vector-valued parameters (permute=True)
    
    sample_dict = fit.extract(permuted=True)
    names = []
    labels = []
    samples = []
    for par, samps in sample_dict.iteritems():
        if par == "lp__":
            loglikes = -samps
        else:
            if len(samps.shape)==1:
                names.append(par)
                labels.append(par)
                samples.append(samps)
            else:
                for i,s in enumerate(samps.T):
                    samples.append(s)
                    names.append(par+str(i+1))
                    labels.append(par+"_"+str(i+1))
    
    return getdist.MCSamples(names=labels, samples=samples, loglikes=loglikes, labels=labels)