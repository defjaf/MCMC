#!/bin/sh

export PYTHONPATH=$HOME/home/proj/stats/MCMC

python driver.py --fdir ./figs_UL2_1 --odir ./out_UL2_1 --idata 0,700 --UL 2 1
python driver.py --fdir ./figs_UL2_2 --odir ./out_UL2_2 --idata 700,1400 --UL 2 1
#python driver.py --fdir ./figs_UL2_3 --odir ./out_UL2_3 --idata 1400,1717 --UL 2 1
