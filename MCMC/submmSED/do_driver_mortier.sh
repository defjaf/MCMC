#!/bin/sh

export PYTHONPATH=$HOME/home/proj/stats/MCMC

export len=49

/usr/local/bin/python driver.py --file ./print_seds_mergers --fdir "./figs_mortier/" --odir "./out_mortier/" --idata $len --UL 0 --format mortier 0
/usr/local/bin/python driver.py --file ./print_seds_mergers --fdir "./figs_mortier/" --odir "./out_mortier/" --idata $len --UL 0 --format mortier 1
/usr/local/bin/python driver.py --file ./print_seds_mergers --fdir "./figs_mortier/" --odir "./out_mortier/" --idata $len --UL 0 --format mortier 2
/usr/local/bin/python driver.py --file ./print_seds_mergers --fdir "./figs_mortier/" --odir "./out_mortier/" --idata $len --UL 0 --format mortier 3

#--idata 175 
#--idata 
#--idata 0,10,1 
#--idata 0,10,1 
#
