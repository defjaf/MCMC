## TODO:

(Some of these are already done???)

* Allow "ancillary" derived parameters/calculations in addition to the MCMC model params (e.g., the ML amplitudes in a marginalized gaussian model -- could store alongside or even as part of MCMC.samples)

* problem when a window is always 0 (or ~0) in a bin for a single experiment? seems to want very large bandpowers? [[iS THIS THE low-ell PROBLEM???]]

* translate/wrap wmap3 likelihood

* polarization in `ClData/` and `binnedCl/` modules

* allow linear combinations of params

* which methods should take the linear combination and which the real parameters? 
  * how about `package`/`unpackage`?


* For sequences, (strings, lists, tuples), use the fact that empty
      sequences are false, so `if not seq` or `if seq` is preferable
      to `if len(seq)` or `if not len(seq)`.
  [ but doesn't work for arrays ]

* just make window fns `win[nbin, ncl, lmax]`

* change mixedCase names to lower\_with\_underscores as per pep?

* get rid of map in favor of list comprehension? (a la python regret)

* get rid of static methods in favor of module-level functions?

* may need to modify `__init__ = blabla` to `def __init__(self, stuff) = self.blabla(stuff)`

* rename `Likelihood` -> `beamLikelihood`?

* `data.get_quadratic_form` method which does A^T N^-1 B or A^T N^-1 A without full inverse if uncorrelated? (Otherwise *likelihood* needs to know which is appropriate)

* allow `MCMC` to take a sequence of likelihoods?
   
* or allow `likelihood` to take a sequence of data?
  * need to wrap `lnLike`, `lnNorm` in something like `[for dat in self.data]`
  * (latter enforces the requirement that they all have the same model, which is crucial, but may be harder to deal with, and requires that every likelihood class do this. or could do with a lnLike1, lnNorm1 which are then wrapped only in the superclass?

* reorganize package/modules



## DONE:


* allow priors from the data? (e.g., reasonable limits on parameter values?)
  * or could subclass the model for specific cases?

* need to make Likelihood take an *instance* of the model class?
  * can't do that -- need to make the xy ranges into *static class* variables.

  * need to link the data, model, likelihood classes more strongly?

  * [done with classmethod]

* how to determine if an object is a vector or an array? `len(asarray(obj))`? 
  * or:
  
              try: l=len(obj)
              except: TypeError: l=1

  * [use `asarray(object).size()` -- is there a better way?]

* generic gaussian data `(d, sigma/sig2/noise_matrix)` subclassed for extra info (x, y for beams, bins for C_l , etc)

----

## `svn2git` notes


Created using [svn2git][] (installed with `sudo gem install svn2git`):

    svn2git svn+ssh://mekon.ph.imperial.ac.uk/Library/Subversion/Repository/MCMC --no-minimize-url --exclude TOIbeams.pickle --authors ../MCMC/authors-transform.txt
    
* Need `--no-minimize` since this repo has several projects at the top level.
* Need `--exclude TOIbeams.pickle` since that file is too big for github.
* Need `--authors ../MCMC/authors-transform.txt` to get full git-style author names 
  * `authors-transform.txt` has been added to root of the project.
  

Push to github at the end (should I use `--all -u` below?):

    git remote add origin https://github.com/defjaf/MCMC.git
    git push -u origin master
    git push --all 
    git push --tags

(Also tried `git push --mirror https://github.com/defjaf/MCMC.git` but this seems to be the official way...)

Note also that `git svn show-ignore` seems to fail, alas.

[svn2git]: https://github.com/nirvdrum/svn2git

----

I previously tried to follow [these instructions][1], but the branch structure seems weird. Could be due to the same structure that required `--no-minimize` above.

    git svn clone svn+ssh://mekon.ph.imperial.ac.uk/Library/Subversion/Repository/MCMC --no-metadata -A authors-transform.txt --stdlayout ~/temp

with the following additions:
* remove large file: `git filter-branch --force --index-filter 'git rm --cached --ignore-unmatch ./TOIbeams.pickle' --prune-empty --tag-name-filter cat -- --all`
 * See <https://help.github.com/articles/removing-sensitive-data-from-a-repository/>
* For some reason, the master branch is called `origin/trunk`, so needed `git branch -m origin/trunk master`
* Push to github at the end: `git push --mirror https://github.com/defjaf/MCMC.git`

[1]: http://john.albin.net/git/convert-subversion-to-git
