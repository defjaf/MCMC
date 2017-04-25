# MCMC

Andrew Jaffe's MCMC engine

Has likelihoods for CMB C_l, nonlinear SED fitting (with normal errors), beam fitting, topology (unfinished?)

## installation

It is probably possible to do this remotely, but easiest is to `git clone` into some directory, 
and then 

    pip install -e .

The `-e` installs in "editable" or "developer" mode which doesn't actually copy the files
into your site-packages, but just leaves them where they are.

It is possible that git will complain about not having `git-lfs` installed. If so, let me know...

There may still be some issues in accessing data files.

The `run` directory shows some recent examples of SED fitting, which is the 
most up-to-date part of the code. There are
other driver files scattered throughout the directory tree which give some other
examples of use, but they may not have been updated for the current directory structure.

## structure of the code

There are separate classes for

* data
* model
* likelihoods

all of which can be subclassed for specific uses.

### TODO:
* Convert to a proper module format -- mostly done
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
