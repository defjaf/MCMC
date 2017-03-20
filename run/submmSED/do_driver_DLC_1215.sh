#!/bin/sh

export PYTHONPATH=${HOME}/home/proj/stats/MCMC

#export hfile=./herus_phot_csv.csv
export hfile=./dat/herus_phot.dat
export lab=DLC_1215_herus
export lin="--lin"

# /Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py ${lin} --format=DLC_2014 --fdir=./figs_${lab} --odir=./out_${lab} --file=./${hfile} 1 1> ${lab}_1.out 2> ${lab}_1.err &

export hfile=./dat/magdis_ulirgs.dat
export lab=DLC_1215_ulirgs
export lin="--lin"

/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py ${lin} --format=DLC_2014 --fdir=./figs_${lab} --odir=./out_${lab} --file=./${hfile} 1 1> ${lab}_1.out 2> ${lab}_1.err &
