"""

Convergence diagnostics

"""

import numpy as np

def gelmanrubin(chains, burn=0, thin=1):
	"""
	Gelman-Rubin convergence diagnostics.
	
	Chains is a (n_chain,) list of MCMC chains.
	
	Each MCMC chain is (n_sample, n_param)  -- need constant n_sample?
	
	""" 
	
	n_chain = len(chains)
	n_samples = [chain.shape[0] for chain in chains]
	n_params = [chain.shape[1] for chain in chains]
	n_param = n_params[0]
	assert all(np == n_param for np in n_params[1:])
	
	print "n_chain = %d" % n_chain
	print "n_param = %d" % n_param
	print "n_samples = ", n_samples
	
	print "Burn-in fraction = %f; thinning = %d" % (burn, thin)
	
	
		
	# var along n_sample (1) axis -- result is (n_chain, n_param)
	# note that the mean is not "debiased" -- divide by (n-ddof)
	within = chains.var(axis=1, ddof=1).mean(axis=0)
	
	## in the ragged-chain formulae (see STAN manual, 54.3)
	within_term = chains.var(axis=1, ddof=0).mean(axis=0)
	between_term = chain_means.var(axis=0, ddof=1) ## nb. divide by (n_chain-1)
	
	totalvar = within_term + between_term
	
	Rhat = np.sqrt(totalvar/within)
	
	return Rhat
	
	
	