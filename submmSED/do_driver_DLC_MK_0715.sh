#!/bin/sh

export PYTHONPATH=${HOME}/home/proj/stats/MCMC

export hfile=./herus_phot_csv.csv
export lab=DLC_MK_2015_2

/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_${lab} --odir=./out_${lab} --file=./${hfile} 0 1> ${lab}_0.out 2> ${lab}_0.err &
/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_${lab} --odir=./out_${lab} --file=./${hfile} 1 1> ${lab}_1.out 2> ${lab}_1.err &
/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_${lab} --odir=./out_${lab} --file=./${hfile} 2 1> ${lab}_2.out 2> ${lab}_2.err &
/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_${lab} --odir=./out_${lab} --file=./${hfile} 3 1> ${lab}_3.out 2> ${lab}_3.err &
/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_${lab} --odir=./out_${lab} --file=./${hfile} 4 1> ${lab}_4.out 2> ${lab}_4.err &
