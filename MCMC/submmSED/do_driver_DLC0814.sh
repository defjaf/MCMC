#!/bin/sh

export PYTHONPATH=${HOME}/home/proj/stats/MCMC

# /usr/local/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2014_2 --odir=./out_DLC_2014_2 --file=./herus_phot.csv 0 1> DLC_2014_0.out 2> DLC_2014_0.err &
# /usr/local/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2014_2 --odir=./out_DLC_2014_2 --file=./herus_phot.csv 1 1> DLC_2014_1.out 2> DLC_2014_1.err &
# /usr/local/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2014_2 --odir=./out_DLC_2014_2 --file=./herus_phot.csv 2 1> DLC_2014_2.out 2> DLC_2014_2.err &
# /usr/local/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2014_2 --odir=./out_DLC_2014_2 --file=./herus_phot.csv 3 1> DLC_2014_3.out 2> DLC_2014_3.err &
# 
/usr/local/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2014_2 --odir=./out_DLC_2014_2 --file=./herus_phot.csv 4 1> DLC_2014_4.out 2> DLC_2014_4.err &
