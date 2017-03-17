#!/bin/sh

export PYTHONPATH=$HOME/home/proj/stats/MCMC

/usr/local/bin/python2.7 driver.py --file=./hrs_planck_sample.dat --format=2 --fdir=DLC0811UL_figs --odir=DLC0811UL_out --idata=65 0
/usr/local/bin/python2.7 driver.py --file=./hrs_planck_sample.dat --format=2 --fdir=DLC0811UL_figs --odir=DLC0811UL_out --idata=65 1
/usr/local/bin/python2.7 driver.py --file=./hrs_planck_sample.dat --format=2 --fdir=DLC0811UL_figs --odir=DLC0811UL_out --idata=65 2
/usr/local/bin/python2.7 driver.py --file=./hrs_planck_sample.dat --format=2 --fdir=DLC0811UL_figs --odir=DLC0811UL_out --idata=65 3

