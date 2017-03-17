#!/bin/sh

export PYTHONPATH=$HOME/home/proj/stats/MCMC

export len=204   #175

/usr/local/bin/python driver.py --file ./pixelfit.dat --fdir "./figs_peel_2/" --odir "./out_peel_2/" --idata $len --UL 0 --format 3 0
/usr/local/bin/python driver.py --file ./pixelfit.dat --fdir "./figs_peel_2/" --odir "./out_peel_2/" --idata $len --UL 0 --format 3 1
/usr/local/bin/python driver.py --file ./pixelfit.dat --fdir "./figs_peel_2/" --odir "./out_peel_2/" --idata $len --UL 0 --format 3 2
/usr/local/bin/python driver.py --file ./pixelfit.dat --fdir "./figs_peel_2/" --odir "./out_peel_2/" --idata $len --UL 0 --format 3 3

#--idata 175 
#--idata 
#--idata 0,10,1 
#--idata 0,10,1 
#
