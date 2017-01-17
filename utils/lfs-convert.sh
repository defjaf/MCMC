#!/bin/bash

### remove *orth.dat from the whole git history (in preparation for use of git-lfs)

### (is there an easier and faster way to do this with wildcards???)

### see https://help.github.com/articles/removing-sensitive-data-from-a-repository/
### and https://help.github.com/articles/configuring-git-large-file-storage/
# 
# for f in likdih/dat/BinaryDihedral2-10orth.dat likdih/dat/BinaryDihedral2-12orth.dat likdih/dat/BinaryDihedral2-14orth.dat likdih/dat/BinaryDihedral2-16orth.dat likdih/dat/BinaryDihedral2-18orth.dat likdih/dat/BinaryDihedral2-20orth.dat likdih/dat/BinaryDihedral2-22orth.dat likdih/dat/BinaryDihedral2-24orth.dat likdih/dat/BinaryDihedral2-26orth.dat likdih/dat/BinaryDihedral2-28orth.dat likdih/dat/BinaryDihedral2-30orth.dat likdih/dat/BinaryDihedral2-32orth.dat likdih/dat/BinaryDihedral2-34orth.dat likdih/dat/BinaryDihedral2-36orth.dat likdih/dat/BinaryDihedral2-38orth.dat likdih/dat/BinaryDihedral2-40orth.dat likdih/dat/BinaryDihedral2-4orth.dat likdih/dat/BinaryDihedral2-6orth.dat likdih/dat/BinaryDihedral2-8orth.dat likdih/dat/BinaryTetrahedral-6orth.dat likoct/dat/BinaryOctahedral-12orth.dat likoct/dat/BinaryOctahedral-16orth.dat likoct/dat/BinaryOctahedral-18orth.dat likoct/dat/BinaryOctahedral-20orth.dat likoct/dat/BinaryOctahedral-24orth.dat likoct/dat/BinaryOctahedral-26orth.dat likoct/dat/BinaryOctahedral-28orth.dat likoct/dat/BinaryOctahedral-30orth.dat likoct/dat/BinaryOctahedral-32orth.dat likoct/dat/BinaryOctahedral-34orth.dat likoct/dat/BinaryOctahedral-36orth.dat likoct/dat/BinaryOctahedral-38orth.dat likoct/dat/BinaryOctahedral-40orth.dat likoct/dat/BinaryOctahedral-8orth.dat liktetr/dat/BinaryTetrahedral-12orth.dat liktetr/dat/BinaryTetrahedral-14orth.dat liktetr/dat/BinaryTetrahedral-16orth.dat liktetr/dat/BinaryTetrahedral-18orth.dat liktetr/dat/BinaryTetrahedral-20orth.dat liktetr/dat/BinaryTetrahedral-22orth.dat liktetr/dat/BinaryTetrahedral-24orth.dat liktetr/dat/BinaryTetrahedral-26orth.dat liktetr/dat/BinaryTetrahedral-28orth.dat liktetr/dat/BinaryTetrahedral-30orth.dat liktetr/dat/BinaryTetrahedral-32orth.dat liktetr/dat/BinaryTetrahedral-34orth.dat liktetr/dat/BinaryTetrahedral-36orth.dat liktetr/dat/BinaryTetrahedral-38orth.dat liktetr/dat/BinaryTetrahedral-40orth.dat liktetr/dat/BinaryTetrahedral-6orth.dat liktetr/dat/BinaryTetrahedral-8orth.dat;
# do
# 	echo ${f}
#     fname=./topology/likelihood/${f}
# 	git filter-branch --force --index-filter \
# 		"git rm --cached --ignore-unmatch $fname" \
# 		--prune-empty --tag-name-filter cat -- --all;
# done

git filter-branch --force --index-filter \
		"git rm --cached --ignore-unmatch ./topology/likelihood/lik\*/dat/\*orth.dat" \
		--prune-empty --tag-name-filter cat -- --all;