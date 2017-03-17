#!/bin/sh

export PYTHONPATH=${HOME}/home/proj/stats/MCMC

/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2015_1 --odir=./out_DLC_2015_1 --file=./herus_phot.csv 0 1> DLC_2015_0.out 2> DLC_2015_0.err &
/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2015_1 --odir=./out_DLC_2015_1 --file=./herus_phot.csv 1 1> DLC_2015_1.out 2> DLC_2015_1.err &
/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2015_1 --odir=./out_DLC_2015_1 --file=./herus_phot.csv 2 1> DLC_2015_2.out 2> DLC_2015_2.err &
/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2015_1 --odir=./out_DLC_2015_1 --file=./herus_phot.csv 3 1> DLC_2015_3.out 2> DLC_2015_3.err &
/Library/Frameworks/Python.framework/Versions/2.7/bin/python driver.py --format=DLC_2014 --fdir=./figs_DLC_2015_1 --odir=./out_DLC_2015_1 --file=./herus_phot.csv 4 1> DLC_2015_4.out 2> DLC_2015_4.err &
