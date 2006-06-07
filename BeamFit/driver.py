import matplotlib
matplotlib.use('Agg')
import pylab

import cPickle

import BeamFit.MAXIPOLBeamData as MP

def main():
    
    dets = [13, 14, 15, 23, 24, 25, 33, 34, 35, 43, 44, 45]
    
    res = {}
    for cols in [(4,5), (6,7)]:
        colname = str(cols[0]).strip()+'_'+str(cols[0]).strip()
        res[cols] = MP.testTOI(nMC=(3000, 100000), useNormalizedBeam=True,
                               noCorrelations=True, fac=None, doBlock=True,
                               cols=None, dets=dets, mapOnly=True, nhits=None,
                               neg=False, rangeScale=None, closeFigs=True,
                               figName='cols_'+colname+'_')
               
    fp = open("all2.pickle", 'wb')  
    
    cPickle.dump(res, fp)   
    
    fp.close()

if __name__ == '__main__':
    main()
