""" 
    read and process 2d beam data (from MAXIPOL) 
    probably need to "duplicate" this module for any significantly different data format and layout
    (e.g., not just different detectors a la MAXIPOL)
"""

from __future__ import division
import sys
import math
import os.path
import pickle
import pylab
import gzip
from itertools import izip, repeat

import numpy
from numpy import  (array, float64, zeros, ones, int32, log, where, exp,
                    arange, asarray, sqrt, minimum, maximum, logical_and)
from numpy import concatenate as cat

from BeamData import BeamData

from ProcessBeamData import setup_sampler, sample1beam, plotter

import getdist

#    sigcut = 0.2 
#    ctscut = 12  
sigcut = 0.02  ## minimum sigma for data -- eliminates weird low-error points
ctscut = 4     ## only keep pixels with this many hits
               ## nb these also apply to TOI data in which case
               ## cts is the number of hits in the 'parent' pixel


def readMAXIPOLdataLuis(filename):
    """ read the data and make an instance of class BeamData """

    ia=[]; ja=[]
    i=[]; j=[]; beam=[]; sig=[]; cts=[]
    for line in open(filename, "r"):
        line = line.strip().split()
        i1, j1, b1, s1, c1 = (int(line[0]), int(line[1]), 
           float(line[2]), float(line[3]), int(line[4]))
        ia.append(i1); ja.append(j1)
        if b1 != 0 and s1 !=0:
            ## only keep pixels with data
            i.append(i1); j.append(j1); beam.append(b1)
            sig.append(s1); cts.append(c1)

    beam = asarray(beam, float64)
    sig = asarray(sig, float64)
    ## map i and j (before deletion) onto (-1,1)
    x = array([2*(ii-min(ia))/(max(ia)-min(ia))-1 for ii in i], float64)
    y = array([2*(jj-min(ja))/(max(ja)-min(ja))-1 for jj in j], float64)

    return BeamData(x, y, beam, sig, cts=cts)
        


def readMAXIPOLdataBrad(filename, day=False, sigcut=0.0, ctscut=0, cols=None,
                        nhits=None, neg=False):
    """ 
    read the data and make an instance of class BeamData 
    use Jeff Collins numbers for Day/Night offset if Day==1

    cols (default (2,3)) gives the column number for data and sigma, respectively
    
    - ignore 'counts' column
    """

    if cols is None: cols=(2,3)

    print "Reading data from columns %d-%d" % tuple(cols)
    
    ngood = 0;
    ncut = 0;
    offsets = {'el': 0.295, 'az': 0.05}   ## degrees (brad)
    #offsets = {'el': 0.289, 'az': 0.065}   ## degrees (Jeff)
    az=[]; el=[]; beam=[]; sig=[]; cts=[]
    if filename.endswith('gz'):
        fileobj = gzip.open(filename, "r");
    else:
        fileobj = open(filename, "r");
    for line in fileobj:
        line = line.strip().split()
#        az1, el1, b1, s1, c1 = (
#           float(line[0]), float(line[1]), float(line[2]), float(line[3]), int(line[4]))
        az1, el1, b1, s1, c1= (
           float(line[0]), float(line[1]), float(line[cols[0]]), float(line[cols[1]]),
           long(line[-1]))
        
        if nhits:
            s1 /= sqrt(c1)

        if day:
            az1 += offsets['az']*60.0   ## arcmin
            el1 += offsets['el']*60.0
            
        if s1>sigcut and c1>ctscut:
            ## only keep pixels with good data
            az.append(az1); el.append(el1); beam.append(b1);
            sig.append(s1); cts.append(c1)
            ngood += 1
        else:
            ncut += 1
            

    fileobj.close()
    
    print 'Data read: ncut=%d, ngood=%d' % (ncut, ngood)

    beam = asarray(beam, float64)
    sig = asarray(sig, float64)
    az = asarray(az, float64)
    el = asarray(el, float64)
    cts = asarray(cts, float64)

    if neg is not False and ((neg is None and beam.mean() < 0) or neg):
        print 'negating data'
        beam = -beam

    return BeamData(az, el, beam, sig, cts=cts)


##MAXI   split from original setup_sampler
def read_data_MAXI(dir=None, files=None, num=None, DayNight=2, LuisBrad=1, cols=(2,3), 
                   sigcut=0.0, ctscut=0, nhits=None, neg=False):
    """
    read the data needed for the sampler for a single beam with data in directory 'dir',
    given by the files in the sequence 'files', or with detector
    'num'.
    
    get the data from column numbers 'cols' from the file
    """
    if dir is None:
        homedir = os.path.expandvars('${HOME}/home')
        if not os.path.exists(homedir):
            homedir = os.path.expandvars('${HOME}')
    
        dir='/'.join((homedir, "cmb/maxipol/beams/"))
    
    if num is not None:
        if LuisBrad == 0: # Luis names
            filen = "binnedb"+str(num).strip()+".txt"
            if DayNight == 0:
                files=["Day/"+filen]
            elif DayNight == 1:
                files=["Night/"+filen]
            elif DayNight == 2:
                files=["Day/"+filen, "Night/"+filen]
        elif LuisBrad == 1: # Brad names
            filen = "brad/b"+str(num).strip()
            suff = 'time_jupiter_maps.txt'
            dnstring = ["day", "night"]
            ## array day is True for daytime, False for nighttime (see read*Brad above)
            if DayNight == 0 or DayNight == 1:
                files=[filen+dnstring[DayNight]+suff]
                day = [not DayNight]
            elif DayNight == 2:
                files=[filen+dnstring[id]+suff for id in [0,1]]
                day = [True, False]
                            
    if files is None:
        files=["Day/binnedb25.txt", "Night/binnedb25.txt"]
        day = [True, False]

    if LuisBrad == 0:
        data = [ readMAXIPOLdataLuis(dir+fil) for fil in files ]
        xyrange = (-1,1)
    elif LuisBrad == 1:
        if num is None and files is not None: 
            day = [False]*len(files)
            print 'Setting day=', day
                
            data = [
                readMAXIPOLdataBrad(dir+fil, d, sigcut=sigcut, ctscut=ctscut, cols=cols, 
                                    nhits=nhits, neg=neg)
                for fil, d in zip(files, day) ]
        xyrange = data[0]   ### set from the zeroth dataset

    return data, xyrange


## this is probably specific to ##MAXI -- just the very inner loop is generic
def sampleall(nruns=2, nMC=(3000, 100000), useNormalizedBeam=True, irun=0,
              noCorrelations=True, fac=None, doBlock=True, nhits=None):
    """
    run the sampler nruns times for the detectors with both Day and Night data
    """
    
    
    
    plotOne = False   ### can't get this to work yet!
    if plotOne:
        pylab.axes()
        pylab.xlim(-100,100)
        pylab.ylim(-100,100)
 
    #dets = [13, 14, 15, 23, 24, 25]   ## LUIS

    if irun==0:
        dets = [12, 13, 14, 15]  # brad, day
        nrow=2; ncol=2
      #  dets = [13, 14] # brad, both
      #  nrow=1; ncol=2
        DayNight=0
    elif irun==1:
        dets = [13, 14, 23, 24, 25, 33, 34, 35, 43, 44, 45 ] # brad, night
        nrow=3; ncol=4
      #  dets = [13, 14] # brad, both
      #  nrow=1; ncol=2
        DayNight = 1
    elif irun==2:    
        dets = [13, 14] # brad, both
        nrow=1; ncol=2
        DayNight=2
    
    reslist = []
    nfig=2
    ntotrun = nruns*nfig
    for run in range(nruns):
        res={}
        for ib, det in enumerate(dets):
            print 'Detector: %d' % det
            fig=pylab.figure(irun*ntotrun+nfig*run)
            if not plotOne:
                ax=fig.add_subplot(nrow, ncol, ib+1)
                ax.cla()

            ## need to explicitly read the data here, now -- how to make generic?
            data, xyrange = read_data_MAXI(num=det, DayNight=DayNight, sigcut=sigcut, ctscut=ctscut)
            like, prop_sigmas, start_params = setup_sampler(data, xyrange,
                                                            useNormalizedBeam=useNormalizedBeam)
            
            res[det] = sample1beam(like, nMC=nMC,  fac=fac,
                                   prop_sigmas=prop_sigmas, start_params=start_params,
                                   noCorrelations=noCorrelations,
                                   doBlock=doBlock)
            if plotOne: 
                pylab.xlim(-100,100)
                pylab.ylim(-100,100)

            fig=pylab.figure(irun*ntotrun+nfig*run+1)
            ax=fig.add_subplot(nrow, ncol, ib+1)
            samples = cat([ s.samples for s in res[det][0] ])
            #samples.transpose()
            ## nb. with numpy, a.transpose doesn't change the array, just gives a new view.
            for var in samples.transpose(): ax.plot(var)
        reslist.append(res)

    return reslist


# probably needs to be split ##MAXI
def testTOI(nMC=(3000, 100000), useNormalizedBeam=True,
            noCorrelations=True, fac=None, doBlock=True, cols=None, dets=None,
            mapOnly=False, nhits=None, neg=None, rangeScale=None, 
            closeFigs=False, figName=None, startCols=None, startParams=None):
    """
    run the sampler  for the detectors with TOI data
    """
    
    reslist = []
    nfig=2
    ntotrun = nfig


 #### AHJ: polz only
 ####if dets is None: dets = [13, 14, 15, 23, 24, 25, 33, 34, 35, 43, 44, 45]
 #####   nrow=4; ncol=3
    
#### needed multipliers for data: -, + are strong requirements, 0 is don't care
#                [33, 34, 43, 45]  [13, 14, 15, 23, 24, 25, 35, 44]
#columns 4-5      -    -   -   -     0  0   0+   0   0  0+   0   0-
#columns 6-7     0-   0-   +   -     -  0    0   +   -   -   0   0
    
    alldets = [33, 34, 43, 45, 13, 14, 15, 23, 24, 25, 35, 44]
    neg45   = [-1, -1, -1, -1, +1, +1, +1, +1, +1, +1, +1, -1]
    neg67   = [-1, -1, +1, -1, -1, +1, +1, +1, -1, -1, +1, +1]

    dfac = {}
    dfac[(2,3)] = dict(izip(alldets, repeat(1)))
    dfac[(4,5)] = dict(izip(alldets, neg45))
    dfac[(6,7)] = dict(izip(alldets, neg67))

    neg1 = neg

    if dets is None: dets = [33, 34, 43, 45]

    nrow = ncol = int(math.sqrt(len(dets)))
    if nrow*ncol < len(dets): ncol += 1
    if nrow*ncol < len(dets): nrow += 1

    pref = 'TOI_polarized/b'
    ident = ['_map','_TOI']
    suff = '.txt.gz'
    if mapOnly: 
        ident = [ident[0]] 

    res={}
    for ib, det in enumerate(dets):

        
        res[det] = []
        startres = []
        filebase = [str(det).strip()+id1 for id1 in ident]

        for fb in filebase:
            fil = pref+fb+suff
            figf = fb
            if isinstance(figName, (str, unicode)):
                figf = figName+figf

            print 'Running: ', fil
            fig=pylab.figure(0)
            ax=fig.add_subplot(nrow, ncol, ib+1)
            ax.cla()
            
            ## this inner part of the loop is probably generic ##MAXI

            try:
                if startCols is not None:
                    if neg is None: neg1 = (dfac[tuple(startCols)][det]<0)

                    data, xyrange = read_data_MAXI(files=[fil],cols=startCols, sigcut=sigcut, ctscut=ctscut,
                                                   nhits=nhits, neg=neg1)
                    like, prop_sigmas, start_params = setup_sampler(
                        data, xyrange,
                        useNormalizedBeam=useNormalizedBeam,
                        rangeScale=rangeScale)

                    startres, ana = sample1beam(like, nMC=nMC, fac=fac, 
                                           prop_sigmas=prop_sigmas,
                                           start_params=start_params,
                                           noCorrelations=noCorrelations,
                                           doBlock=doBlock)
                                           
                    mod = like.model
                    orig_params = mod.package(ana[0])     ## startres[-1].mean()
                    orig_sigmas = mod.package(3*ana[1])   ## startres[-1].stdev()

                    print ("Start: " + mod.fmtstring) % tuple(mod.unpackage(orig_params))
                    print ("Sigma: " + mod.fmtstring) % tuple(mod.unpackage(orig_sigmas))
                    
                ## need to run this to get the correct likelihood.
                ##   therefore may need to adjust prior ranges

                if neg is None: neg1 = (dfac[tuple(cols)][det]<0)
                data, xyrange = read_data_MAXI(files=[fil],cols=cols, sigcut=sigcut, ctscut=ctscut,
                                               nhits=nhits, neg=neg1)
                like, prop_sigmas, start_params = setup_sampler(
                    data, xyrange,
                    useNormalizedBeam=useNormalizedBeam,
                    rangeScale=rangeScale)
                
                if startCols is not None:
                    prop_sigmas = orig_sigmas
                    start_params = orig_params
                    
                if startParams is not None:
                    start_params=startParams
                
                #prev_fig = pylab.gcf()    
                #pylab.figure(4)
                #ll, xx, yy = get_likelihood_grid(like, start_params) 
                ##return ll,xx,yy
                #pylab.figure(prev_fig)

                res[det].append(sample1beam(like, nMC=nMC, prop_sigmas=prop_sigmas,
                                       start_params=start_params, fac=fac, 
                                       noCorrelations=noCorrelations, doBlock=doBlock))
                
                if figName:
                    fig.savefig(figf+str(fig.number).strip()+'.png')

                sys.stdout.flush()
                fig=pylab.figure(1)
                ax=fig.add_subplot(nrow, ncol, ib+1)
                samples = cat([ s.samples for s in res[det][-1][0] ])

                for var in samples.transpose(): ax.plot(var)
                if figName:
                    fig.savefig(figf+str(fig.number).strip()+'.png')

                fig=pylab.figure(2)
                getdist.histgrid(res[det][-1][0][-1])

                if figName:
                    fig.savefig(figf+str(fig.number).strip()+'.png')

            except None:
                print "Unexpected error:", sys.exc_info()[0]
                print "... when running ", fil, fb

            if closeFigs: pylab.close('all')

    pylab.show()
    return res


def saveres(reslist, file=None):
    """
    pickle the results; can't save the actual sampler class since
    you can't pickle all its members
    
    NO LONGER NEEDED; can just pickle the MCMC instances themselves
    """
    newres=[]
    for resrun in reslist:
        for det, res in resrun.iteritems():
            resrun[det]=(res[0][-1].samples, res[1])  ## save the 'ana' element & samples
            ## should save in a form closer to the original???
        newres.append(resrun)
    if file is not None:
        fp = open(file, 'w')
        pickle.dump(newres, fp)
        fp.close()
    return newres


def makereport(reslist, file=sys.stdout, hasTOI=False, hasRuns=False):
    """
    print out a 'report' of the MCMC run
    """

    if not hasRuns: reslist=[reslist]
    
    for irun, resrun in enumerate(reslist):
        file.write("Run: %d\n" % irun)
        for det, res in resrun.iteritems():
            file.write("%d" % det)

            try:
                if not hasTOI: 
                    val = res[1][0]
                    sig = res[1][1]
                else:
                    val = res[-1][1][0]
                    sig = res[-1][1][1]
                
                for (v, s) in zip(val, sig):
                    file.write("   %f +- %f" % (v, s))
            except:
                print "\n... when running ", irun, det,
                print "Unexpected error:", sys.exc_info()[0]
                    
            file.write("\n")

    
