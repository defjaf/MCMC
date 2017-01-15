# MCMC

Andrew Jaffe's MCMC engine

Has likelihoods for CMB C_l, nonlinear SED fitting (with normal errors), beam fitting, topology (unfinished?)

----

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

[svn2git]: https://github.com/nirvdrum/svn2git
