# MCMC

Andrew Jaffe's MCMC engine

Has likelihoods for CMB C_l, nonlinear SED fitting (with normal errors), beam fitting, topology (unfinished?)


### TODO:
* Convert to a proper module format
  * How/where to store data and outputs in module?
  * Move up to "above" module dir?
  * probably need to be careful of relative imports
* Split off `svn2git` stuff from below


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
