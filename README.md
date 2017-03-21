# MCMC

Andrew Jaffe's MCMC engine

Has likelihoods for CMB C_l, nonlinear SED fitting (with normal errors), beam fitting, topology (unfinished?)


### TODO:
* Convert to a proper module format
  * How/where to store data and outputs in module?
    * done in setup.py
    * MANIFEST.in lists non-python files for module installation 
      * not needed
    * nb. differentiate between what's on github and what gets installed
    * Move up to "above" module dir?
      * see new "run" directory
  * need to be careful of relative imports?
    * Done (mostly?)
  * move shell driver code and outputs (and data??) out of the package into "run" or "work" dir?
    * started (for submmSED)
      * currently symlink to MCMC/submmSED/dat directory
  * move stan code and runs to a different package or submodule?
