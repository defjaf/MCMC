
import pickle
import getdist_ahj
import pylab
import numpy

suff = ''  # '_1'

def getpickle(file=None):
    if file is None:
        file = 'topology/wmap06res'+suff+'.pickle'
    fp = open(file)
    return pickle.load(fp)

def main(file = None):

    res = getpickle(file=file)

    for i, (topo, mcmc) in enumerate(res.items()):
        pylab.figure(i)
        print('topology: %s' % (topo))
        pylab.figtext(0.75, 0.75, topo)
        
        getdist_ahj.histgrid(mcmc[0][-1].samples, [0,1,2,3,4])
        if i>0: break
        #pylab.savefig(topo+suff+'.png')
        