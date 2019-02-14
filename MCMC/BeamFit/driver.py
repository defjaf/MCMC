import matplotlib
#matplotlib.use('TkAgg')
import pylab

import pickle

from . import MAXIPOLBeamData as MP

def main_TOI():
    
    dets = [13, 14, 15, 23, 24, 25, 33, 34, 35, 43, 44, 45]

    #cols_list = [(2,3)]
    cols_list = [(4,5), (6,7)]

    mapOnly = True 
    res = {}
    for cols in cols_list:
        colname = str(cols[0]).strip()+str(cols[1]).strip()
        res[cols] = MP.testTOI(nMC=(3000, 100000), useNormalizedBeam=True,
                               noCorrelations=True, fac=None, doBlock=True,
                               cols=cols, dets=dets, mapOnly=mapOnly, nhits=None,
                               neg=False, rangeScale=None, 
                               closeFigs=True,
                               figName='QUBeams.out/cols_'+colname+'_')
               
    fp = open("all2.pickle", 'wb')  
    
    pickle.dump(res, fp)   
    
    fp.close()

def main_maps():

    dets = [13, 14, 15, 23, 24, 25, 33, 34, 35, 43, 44, 45]

    cols_list = [(2,3)]
    #cols_list = [(4,5), (6,7)]:

    mapOnly = True 
    res = {}
    for cols in cols_list:
        colname = str(cols[0]).strip()+str(cols[1]).strip()
        res[cols] = MP.testTOI(nMC=(1000, 10000), useNormalizedBeam=True,
                               noCorrelations=True, fac=None, doBlock=True,
                               cols=cols, dets=dets, mapOnly=mapOnly, nhits=None,
                               neg=False, rangeScale=None)

    fp = open("allmaps.pickle", 'wb')  

    pickle.dump(res, fp)   

    fp.close()

if __name__ == '__main__':
    main()
