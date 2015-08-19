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
    
    