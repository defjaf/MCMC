"""

Convergence diagnostics

"""

import numpy as np

def gelmanrubin(chains_complete, burn=0, thin=1, verbose=False):
    """
    Gelman-Rubin convergence diagnostics.
    
    chains_complete is a (n_chain,) list of MCMC chains.
    
    Each MCMC chain is (n_sample, n_param)  -- need constant n_sample?
    
    in the analysis, remove the first burn samples, and use every stride samples
    
    """ 
    
    n_chain = len(chains_complete)
    n_samples = [chain.shape[0] for chain in chains_complete]
    n_params = [chain.shape[1] for chain in chains_complete]
    n_param = n_params[0]
    assert all(np == n_param for np in n_params[1:])
    
    if verbose:
        print "n_chain = %d" % n_chain
        print "n_param = %d" % n_param
        print "n_samples = ", n_samples
    
        print "Burn-in fraction = %f; thinning = %d" % (burn, thin)
    
    chains = [c[burn::stride] for c, ns in zip(chains_complete, n_samples)]
        
    ### can't index by all three indices since it's a ragged list of arrays
    ### use list comprehensions for the first var/mean operations below
            
    # var along n_sample axis via list -- result is (n_chain, n_param)
    # followed by mean along chain axis -- result is (nparam,)
    # note that var is "debiased" -- divide by (n-ddof)
    within = np.array([np.var(c, axis=0, ddof=1) for c in chains]).mean(axis=0)
    
    ## in the ragged-chain formulae (see STAN manual, 54.3)
    ## the within term isn't debiased, unlike above.
    within_term = np.array([np.var(c, axis=0, ddof=0) for c in chains]).mean(axis=0)
    between_term = np.array([np.mean(c, axis=0) for c in chains]).var(axis=0, ddof=1) ## nb. divide by (n_chain-1)
    
    totalvar = within_term + between_term
    
    Rhat = np.sqrt(totalvar/within)
        
    return Rhat


#### need to make a version which runs on output from my sampler

def gelmanrubin_MCMC(MCMCs, burn=0, stride=1, verbose=False):
    """
    Gelman-Rubin convergence diagnostics.
    
    MCMCs is a (n_chain,) list of MCMC objects.
        
    in the analysis, remove the first burn*nsample samples, and use every thin samples

    *** use the saved/computed parameter means and variances    
    
    TODO: also evaluate for the chi^2 or lnPr at max?
    
    """ 
    chain_means = np.array([m.mean(burn, stride) for m in MCMCs])
    chain_vars = np.array([np.power(m.stdev(burn, stride), 2) for m in MCMCs])
    
    #### the actual ragged-chain definition uses ddof=1 for the "within-chain" calculation, 
    ####    but ddof=0 for the within-chain part of the total variance. 
    ####    Ignore that here, since it's much easier to code this way. (Nsamp>>1)
    within = chain_vars.mean(axis=0)            ### the mean of the variances
    between = chain_means.var(axis=0, ddof=1)   ### the variance of the means

    ### there can be zero-variance (fixed) parameters
    idxs = within>0
    within = within[idxs]
    between = between[idxs]

    Rhat = np.sqrt((within+between)/within)

    return Rhat
