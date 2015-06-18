#!/usr/bin/env python
# encoding: utf-8
"""
driver.py

Created by Andrew H. Jaffe on 2010-05-22.
Copyright (c) 2010 Imperial College London. All rights reserved.
"""

from __future__ import division

import sys
import os

import getopt

import operator
import cPickle as pickle

import numpy as np

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("AGG")
    plottype = "png"
    
    #print "NOT USING LaTeX"; matplotlib.rc('text', usetex=False)  ## memory leak with latex????
    
import matplotlib.pyplot as plt

import MCMC
import likelihood
import data
import model
import M31model
import getdist_ahj as getdist

# import joblib


#### pattern after test_binnedCl.py and/or BeamFit/driver.py BeamFit/MAXIPOLBeamData.py

fname_DLC = "./submmSED.txt"
fname_ERCSC = "./submmSED/ercsc_iifscz.txt"
fname_MRR_old = "./submmSED/ERCSCalliifscz4550850.dat"
fname_MRR = "./submmSED/ERCSCiifsczbg.dat"
fname_Peel = "./submmSED/M31/pixelfit.dat"
fname_Mortier = "./submmSED/print_seds_mergers"
fname_M31 = "./submmSED/M31/M31Flux-v2.dat"
fname_DLC_2014 = "./submmSED/herus_phot.csv"
delnu = 1763



wavelength = True ### Planck format

### TODO: add calculation of the likelihood/posterior of the posterior mean params
def main(filename=fname_MRR, i=None, rotateParams=False, onecomponent=True, getNorm=True, start=None, sigmas=None, 
         nMC=(10000,10000), nDerived=None, noPlots=False, fig0=0, savefig=False, retMCMC=True,
         opticallyThick=False,
         random=True, randomrestart=False,
         fdir = "./", logplot=True, DLC_ul=False, check=None, next0=True, format=0, filetype='MRR'):
        
        
    speedOfLight = 299792.  ## micron GHz
    nu2, nu1 = speedOfLight/8.0, speedOfLight/1000.0 ### rest-frame microns -> GHz
    
    noHist = False  
    
    ret = []
        
#    start = np.array(start) + 0.5*np.array(sigmas)*np.random.randn(len(start))
#    print "randomizing start:", start
            
    ## read the data
    if filetype.upper() == 'DLC':
        alldata = data.readfluxes_DLC(filename, format=format, delnu=delnu)   ### nb. format=3 is PEEL
    elif filetype.upper() == 'MRR':
        if DLC_ul==1:
            print "Removing 12micron and 217GHz; 25micron UL"
            alldata = data.readfluxes_MRR(filename, IRAS_ignore=[0], Planck_ignore=[3], DLC_ul=True, next0=next0)
        elif DLC_ul==2:
            print "Removing 12micron, 353 GHz, 217GHz; 25micron UL"
            alldata = data.readfluxes_MRR(filename, IRAS_ignore=[0], Planck_ignore=[2,3], DLC_ul=True, next0=next0)
        else:
            alldata = data.readfluxes_MRR(filename)
    elif filetype.upper() == 'PEEL':    ### can also read with DLC format=3
        alldata = data.readfluxes_peel(filename, delnu=delnu)
    elif filetype.upper() == 'MORTIER':
        alldata = data.readfluxes_mortier(filename)
    elif filetype.upper() == "M31":
        alldata = data.readfluxes_M31()
    elif filetype.upper() == "DLC_2014":
        alldata = data.readfluxes_DLC_2014()
    else:
        alldata = data.readfluxes_ERCSC_TopCat(filename)
        
    if i is None:
        idata = range(len(alldata))
    elif not operator.isSequenceType(i):
        idata = [i]
    else:
        idata = i
        
    for i in idata:

        try:
            dat = alldata[i]
            name = dat.name + "(z=%5.3f)" % (dat.z)
            singleObject = True
        except TypeError:
            dat = [alldata[j] for j in i]
            name = " + ".join(str(int(d.name)) for d in dat)
            singleObject = False
        except IndexError:
            print "Got IndexError -- no more objects at i=%d" % i
            break   ## or continue?
        
        print "Object[s] %s" % name
    
        ## initialize the model (with a class, not an object)
        
        if filetype.upper() == "M31":
            mod = M31model.M31model
            like = likelihood.SEDLikelihood_normalized(data=dat, model=mod)
        elif opticallyThick:
            mod = model.submmModel1_opticallythick_logA
            like = likelihood.SEDLikelihood_normalized(data=dat, model=mod)   
        elif getNorm:
            if onecomponent:
                mod = model.submmModel1_normalized_logA
            else:
                mod = model.submmModel2_normalized_logA
            like = likelihood.SEDLikelihood_normalized(data=dat, model=mod)        
        else:
            if onecomponent:
                mod = model.submmModel1
                like = likelihood.SEDLikelihood1(data=dat, model=mod)
            else:
                mod = model.submmModel2            
                like = likelihood.SEDLikelihood2(data=dat, model=mod)
        
        if start is None: 
            start_params = np.asarray(mod.startfrom(random=random))
        else:
            start_params = np.asarray(start)
        if sigmas is None:
            prop_sigmas = start_params/4.0
        else:
            prop_sigmas = sigmas
        
        if nDerived is not None:
            like.nDerived = nDerived
    
        mcmc, ana = MCMC.sampler(like, nMC, prop_sigmas, start_params, plotter=None,
                            fac=None, noCorrelations=True, doBlock=True, rotateParams=rotateParams,
                            randomrestart=randomrestart)

        if not noPlots:

            fig = plt.figure(fig0)
            lnLike = []
            plt.subplots_adjust(wspace=0.3, hspace=0.25)
            maxlnLike, maxLikeParams = getdist.histgrid(
                                         mcmc[-1], noPlot=noHist, params=np.where(np.array(sigmas)>0)[0],
                                         burn=0.2, stride=1)
            plt.subplots_adjust()

            #plt.suptitle(name)
    
            if savefig:
                try:
                    fname = fdir+name + savefig
                except TypeError:
                    fname = fdir+name
                    
                try:     ### move close to finally clause? 
                    fig.savefig(fname+"_0."+plottype)
                    plt.close(fig)
                except Exception as err:  ## dangerous -- catch anything!
                    print "CAUGHT Exception!! -- WHILE SAVING %s" % fname+"_0.png"
                    print "Error Info: "
                    print type(err)     # the exception instance
                    print err.args      # arguments stored in .args
                    print err
    
            fig=plt.figure(fig0+1)
            params = ana[0]
            meanmod = mod(*params)
            meanlnProb = like.lnLike(params)
            mean_chi2 = like.chi2(params)
            print "ln Pr of mean = %f" % meanlnProb
            print "chi2(mean) = %f"% mean_chi2

            MLmod = mod(*maxLikeParams)
            try:
                meanmod.plot(dat, wavelength=wavelength, logplot=logplot)
                MLmod.plot(dat, wavelength=wavelength, logplot=logplot)
                # plt.suptitle(name)    
            except AttributeError:
                pass
            if savefig:
                try:     ### move close to finally clause? 
                    fig.savefig(fname+"_1."+plottype)
                    plt.close(fig)
                except Exception as err:   ## dangerous -- catch anything!
                    print "CAUGHT Exception!! -- WHILE SAVING %s" % fname+"_1.png"
                    print "Error Info: "
                    print type(err)     # the exception instance
                    print err.args      # arguments stored in .args
                    print err
                
                
            fig0 += 2
    
        ret_i= (ana, (maxlnLike, maxLikeParams, meanlnProb), name) ### meanlnProb added

        if singleObject:
            ### collect further information to return
            ret_i += (dat.z,)
            ret_i += (zip(dat.d, dat.sig),)
            ret_i += (MLmod.flux(nu1, nu2),)
    
    
        if retMCMC:
            ret.append((mcmc,)+ret_i) 
        else:
            del mcmc
            ret.append(ret_i)
            
        if check is not None:
            with open(check, 'w') as f:
                pickle.dump(ret, f)
            
            
    
    if len(ret) == 1:
        ret = ret[0]
        
    return ret
    
    
def recover():
    """ to recover from checkpoint files -- UNFINISHED, but note that postprocess now works on checkpoints"""
    with open("./submmSED/out_next0/check0.npy") as f:
        ret0 = pickle.load(f)
    with open("./submmSED/out_next0/check1.npy") as f:
        ret1 = pickle.load(f)
        
    

# idata =[i*50 for i in range(6,29)]
# nMC = (15000,100000)

#idata =[i*25 for i in range(12,57)]
nMC = (15000,300000)
#    fil = "./ercsc_iifscz.txt"
idata = range(0,1717) #None   #[0,300,700] #
# fil = "./ERCSCalliifscz4550850.dat"
fil = "./ERCSCiifsczbg.dat"


def many(which = range(4), idata=idata, nMC = nMC, fil=fil, fdir="./", cdir="./", next0=True, **keywords):

    print "Using file %s" % fil
    
    print "keywords:", keywords
    
    ret1 = ret2 = ret3 =  ret4 = []

    ### proposition sigmas ###
    sA, sB, sT, snu = 4, 8, 8, 10   ## was 1,2,2

    if 0 in which:
        print "Two-Component beta = 2"
        ret1 = main(fil, getNorm=True, i = idata, 
                    start=(1,2.,10,0.1,2.,20), sigmas=(sA,0,sT)*2, retMCMC=False,
                    nMC=nMC, onecomponent=False, fig0=0, savefig="_2comp_b2", fdir=fdir,
                    check=cdir+"check0.npy", next0=next0, **keywords)

    if 1 in which:
        print "One-Component"
        ret2 = main(fil, getNorm=True, i = idata, 
                    start=(1,2.,10), sigmas=(sA,sB,sT), retMCMC=False,
                    nMC=nMC, onecomponent=True, fig0=100, savefig="_1comp", fdir=fdir,
                    check=cdir+"check1.npy", next0=next0, **keywords)
                
    if 2 in which:
        print "One-Component beta = 2"
        ret3 = main(fil, getNorm=True, i = idata, 
                    start=(1,2.,10), sigmas=(sA,0,sT), retMCMC=False,
                    nMC=nMC, onecomponent=True, fig0=200, savefig="_1comp_b2", fdir=fdir,
                    check=cdir+"check2.npy", next0=next0, **keywords)
                
                
    if 3 in which:
        print "Two-Component"
        ret4 = main(fil, getNorm=True, i = idata, 
                    start=(1,2.,10,0.1,2.,20), sigmas=(sA,sB,sT)*2, retMCMC=False,
                    nMC=nMC, onecomponent=False, fig0=0, savefig="_2comp", fdir=fdir,
                    check=cdir+"check3.npy", next0=next0, **keywords)

    if 4 in which:
        print "Optically thick"
        ret4 = main(fil, getNorm=True, i = idata, 
                    start=(1,2.,10, 100.0), sigmas=(sA,sB,sT,snu), retMCMC=False,
                    nMC=nMC, onecomponent=False, opticallyThick=True, fig0=0, savefig="_thick", fdir=fdir,
                    check=cdir+"check4.npy", next0=next0, **keywords)


    return ret1, ret2, ret3, ret4
    

def postprocess(dirname="./", multiple=None, check=False, nodat=False):
    """
    process the pickle files returned by main() and many() routines.
    Can also deal with the slightly-different format of the checkpoint files.
    """    
    #### TODO: save more information for DLC (z, fluxes, error, evidence, total fluxes) DONE
    #### allow combining different pickle files DONE
    
    #### TODO: fix the normalization on the evidence calculations (priors are currently incorrect)
    ####        can use Savage-Dickey?
    
    
    nrun = 5
    
    ret = [[] for i in range(nrun)]
    idxs = [[0,2,3,5], [0,1,2], [0,2], [0,1,2,3,4,5], [0,1,2,3]]  ## penultimate one is for full 2T-2beta fit
    nt = [2,1,1,2,1]  ## number of temperature components
    
    
    if not multiple:
        dirname = [dirname]
    
    for i in range(nrun):
        
        for dirn in dirname:
        
            if check:
                fname = dirn+"check%d.npy" % i                        
            else:
                fname = dirn+"out_[%d].pickle" % i        
            print fname
            try:
                with open(fname) as f:
                    ret0 = pickle.load(f)
            except IOError:
                continue
            
            ix = idxs[i]
            if not check: 
                ret0 = ret0[i]
            nobj = len(ret0)
            npar = len(ix)
            try:
                ndat = len(ret0[0][4])
            except IndexError:
                ndat = 0
            if nodat:
                ndat = 0
                
            print 'nobj, npar, ndat = ', nobj, npar, ndat
                            
            dt = np.dtype([
                ('name', 'S21'),
                ('mean', np.float, (npar,)), 
                ('sig', np.float, (npar,)), 
                ('covar', np.float, (npar,npar)), 
                ('ML', np.float),
                ('ev', np.float),
                ('MLpar', np.float, (npar,)),
                ('MeanL', np.float),
                ('evMean', np.float),
                ('dlnLike', np.float),
                ('z', np.float),
                ('dat', np.float, (ndat,2)),
                ('flux', np.float, (nt[i],))
            ])
        
            ret_i = np.empty(nobj, dtype = dt)
        
            ### each of ret[i] has format
            ## [ <for each object>... 
            ##   (
            ##     ([parameter means], [parameter sigmas], [normalized covar]),
            ##     (ML, [ML params], meanL),  !!! meanL added recently
            ##     name    !!! added in recent versions
            #    )
            ## ...]
 
            for iobj, obj in enumerate(ret0):
                ret_i[iobj]['name'] = obj[2]
                ret_i[iobj]['mean'][:] = np.array(obj[0][0])[ix]
                ret_i[iobj]['sig'][:] = np.array(obj[0][1])[ix]
                covar = np.array(obj[0][2])[:,ix][ix]
                ML = obj[1][0]
                ret_i[iobj]['covar'][:,:] = covar
                ret_i[iobj]['ML'] = ML
                ret_i[iobj]['ev'] = ML + 0.5*np.linalg.det(covar) + npar*0.5*np.log(2*np.pi)
                ret_i[iobj]['MLpar'][:] = obj[1][1][ix]
                try:
                    meanL = obj[1][2]
                    ret_i[iobj]['MeanL'] = meanL
                    ret_i[iobj]['evMean'] = meanL + 0.5*np.linalg.det(covar) + npar*0.5*np.log(2*np.pi)
                    ret_i[iobj]['dlnLike'] = ML-meanL
                    
                    ### new DLC information
                    ret_i[iobj]['z'] = obj[3]
                    if ndat>0:
                        ret_i[iobj]['dat'][:,:] = np.array(obj[4])[:,:]
                    ret_i[iobj]['flux'][:] = np.array(obj[5])[:]  ### AHJ: PROBLEM HERE-- fixed with [:]
                except IndexError:
                    pass
                    
                    
            ret[i].append(ret_i)
            
        # concatenate different output files (in different directories as written)
        if ret[i]:
            ret[i] = np.concatenate(ret[i])
    
    return ret


def writeTabAll(ret123, fbase, ext='.npy', dirname=None, check=False, nodat=False):
    if dirname is not None:
        ret123 = postprocess(dirname, check=check, nodat=nodat)
        
    for i, r in enumerate(ret123):
        if len(r)>0:
            fname = fbase + str(i) + ext
            writeTab(r, fname, nodat=nodat)
        
    return ret123
        

def writeTab(ret, fname, names=None, nodat=False):
    """ write the output of the postprocess function to a text file 
        run separately on each of the elements of postprocess()
        if nodat, don't include the data (needed for datasets with different number of data per object
     """
     #### TODO: write more information for DLC (z, fluxes, error, evidence, total fluxes)
     
    try:
        anames = ret['name']
    except ValueError:
        anames = np.array(names)

    nn = ret.shape[0]
    npar = len(ret['MLpar'][0])
    ndat = len(ret['dat'][0]) if not nodat else 0
    
#    try:
    nt = ret['flux'].shape[1]
    
    alls = np.hstack([anames.reshape(nn,1), ret['z'].reshape(nn,1),
          ret['MLpar'], ret['mean'], ret['sig'], ret['dlnLike'].reshape(nn,1), 
          ret['ev'].reshape(nn,1), ret['evMean'].reshape(nn,1),
          ret['dat'].reshape(nn,-1), ret['flux'].reshape(nn,-1)    ### remove these two lines for old files
          ])
   # except ValueError:
    #      alls = np.hstack([anames.reshape(nn,1), ret['z'].reshape(nn,1),
    #            ret['MLpar'], ret['mean'], ret['sig'], ret['dlnLike'].reshape(nn,1), 
    #            ret['ev'].reshape(nn,1), ret['evMean'].reshape(nn,1)
    #            ])
    #      nt = 0
         
    print 'nn, npar, ndat, nt = ', nn, npar, ndat, nt
        
    hdr = ['Name', 'z']
    for i in range(npar):
        hdr.append("ML param %d" % i)
    for i in range(npar):
        hdr.append("Mean param %d" % i)
    for i in range(npar):
        hdr.append("sigma param %d" % i)
    hdr.append("dlnLike")
    hdr.append("evidence1")
    hdr.append("evidence2")
    if not nodat:
        for i in range(ndat):
            hdr.append("flux %d" % i)
            hdr.append("sigflux %d" % i)
    for i in range(nt):   ### at least one of these is printed incorrectly!!!
        hdr.append("greybody flux %d" % i)
        
    nhead = len(hdr)
    hdr = ("%21s "*nhead) % tuple(hdr)
    with open(fname, 'w') as f:
        f.write(hdr + '\n')
        np.savetxt(f, alls, fmt="%21s", delimiter=' ')
        
        
def plotter(sampler):
    """
    make a plot from the results of the sampler
    (could probably make this more generic if plotMod is a member fn of the model)
    for now, only plot the *first* dataset (since for this case that's the 'good' data)
    
    TODO: calculate stats for the amplitude here?
    
    NOT USED IN ABOVE YET
    """

    mod = sampler.like.model
    data = sampler.like.data

    params = where(sampler.prop.sigmas>0)[0]


    ntot = len(sampler.samples)
    stride = 1 #max(1,ntot//sampler.naccept)
    ana = MCMC.chain_analyze(sampler.samples[(ntot//5)::stride,:], params=params)
    vals = sampler.like.model.package(ana[0])
    sigs = sampler.like.model.package(ana[1])

    ### need to take shape into account
    vals = mod.bandpowers(vals)
    sigs = mod.bandpowers(sigs)


    #vals = vals*mod.ellctr*(mod.ellctr+1)/(2*math.pi)
    #sigs = sigs*mod.ellctr*(mod.ellctr+1)/(2*math.pi)


    print vals  #sampler.like.model.fmtstring % tuple(ana[0])
    print sigs  #sampler.like.model.fmtstring % tuple(ana[1])
    print ana[2]


    #m = mod(vals); c1 = m()[0]; ell = N.arange(c1.size)
    #c = c1*ell*(ell+1)/(2*N.pi)
    #pylab.plot(ell, c)
    
# def simul():
#     """ supposed to do embarassingly parallel python launches, but doesn't appear to work... """
#     par = joblib.Parallel(n_jobs=3, verbose=1)
#     ret123 = par(joblib.delayed(many([w])) for w in range(3))
#     with open("out_123.pickle", 'w') as f:
#         pickle.dump(ret123, f)
#         
#     return ret123


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def mainmain(argv=None):
    """
    run submmSED MCMC
    
    python driver [options] which0 [which1 ...]
    --file: input data filename
    --fdir, -f: figure directory
    --odir, -o: output numpy pickle files
    --idata, -i: comma/space-separated list in python slice format (start, stop, step) of which data to use
    --UL: use DLC upper-limit calculation
    --format: DLC file format (default 0; see data.py)
    which = 0: 2 comp b=2
            1: 1 comp floating b
            2: 1 comp b=2
            3: 2 comp floating b
    """
    
    
    randomrestart = True
    random = True
    
    # fdir = "./figs_MRR_UL_XX/"
    # odir = "./out_MRR_UL_XX/"
    fdir = "./figs_DLC2_1/"
    odir = "./out_DLC2_1/"
    
    fdir = "./figs_Planckfinal/"
    odir = "./out_Planckfinal/"
    DLC_ul = True
    filetype = 'MRR'
    which = []
    datslice=None
    format=None
    filename = fil
    next0 = False ### usually want True, but need this for indexing by the line numbers in the file
    
    longopts = ["help", "fdir=", "odir=", "idata=", "UL=", "format=", "file="]
    shortopts = "hf:o:i:"
    if argv is None:
        argv = sys.argv

    try:
        try:
            opts, args = getopt.getopt(argv[1:], shortopts, longopts)
        except getopt.GetoptError as msg:
             raise Usage(msg)

        print "normal behaviour"
        for o, a in opts:
            if o in ("-h", "--help"):
                print mainmain.__doc__
                sys.exit(0)
            elif o in ("--fdir", "-f"):
                fdir = a
            elif o in ("--odir", "-o"):
                odir = a
            elif o in ("-i", "--idata"):
                datslice = [int(i) for i in a.replace(","," ").split()]
                idata = range(*datslice)
            elif o in ["--UL"]:
                DLC_ul = int(a)
            elif o in ["--format"]:
                try:
                    format = int(a)
                except ValueError:
                    format = None
                    filetype = a
            elif o in ["--file"]:
                filename = a
                
                
            if not fdir.endswith('/'): fdir+='/'
            if not odir.endswith('/'): odir+='/'

            try:
                os.mkdir(odir)
            except OSError:
                pass
            try:
                os.mkdir(fdir)
            except OSError:
                pass
                
        if format is not None:
            filetype = 'DLC'
        if format==4:
            filetype = 'DLC_2014'
            DLC_UL = 0
            format = 0
        else:
            format = 0 ### placeholder

        if datslice is None:
            if filetype == 'DLC':
                idata = range(1717)
            elif filetype == 'DLC_2014':
                idata = range(43)
        
        print "filename: %s" % filename 
        print "fig dir: %s" % fdir
        print "out dir: %s" % odir
        print "data range", datslice
        print "DLC_UL = %d" % DLC_ul
        print "format = %d" % format
        print "filetype = %s" % filetype
                
        # process arguments
        for s in reversed(args):
            try:
                which.append(int(s))
            except ValueError:
                break
        print "which=", which
        ret = many(which, fdir=fdir, DLC_ul=DLC_ul, filetype=filetype, cdir=odir, idata=idata, next0=next0, 
                   fil=filename, format=format, random=random, randomrestart=randomrestart)
        with open(odir+"out_"+"".join(str(which).split(' '))+".pickle", 'w') as f:
            pickle.dump(ret, f)
            
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(mainmain())
# 
# 
# ### TODO: allow setting directory and idata from the command line
# if __name__ == '__main__':
#     # fdir = "./figs_MRR_UL_XX/"
#     # odir = "./out_MRR_UL_XX/"
#     fdir = "./figs_XX3/"
#     odir = "./out_XX3/"
#     DLC_ul = True
#     which = []
#     for s in reversed(sys.argv):
#         try:
#             which.append(int(s))
#         except ValueError:
#             break
#     print "which=", which
#     ret = many(which, fdir=fdir, DLC_ul=DLC_ul, cdir=odir)
#     with open(odir+"out_"+"".join(str(which).split(' '))+".pickle", 'w') as f:
#         pickle.dump(ret, f)
    
    
    

