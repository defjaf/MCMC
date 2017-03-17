The IIFSCz Catalogue has been revised to correct the fluxes for extended sources,
and to make some improvements to the infrared template fitting.
A further revision has been made to improve the accuracy of the luminosities
at low redshift (24/2/10).
The adopted position for 6DF associations has been revised to be the 6DF position
(rather than the IRAS position, as previously). (18/6/10)

The subset of hyperluminous galaxies (L_ir>10^13 Lo) is IIFSCz_hyper.dat (with same format).
(16/4/10)

The columns of IIFSCz catalogue (IIFSCz_v5.dat.gz, 60303 sources) are:

nameIRAS, ra, dec, posFlag,
S12, S25, S60, S100, nq1, nq2, nq3, nq4,
am1, am2, am3, am4, am5, am6, am7, am8, em1, em2, em3, em4, em5, em6, em7,
em8, photFlag,
FINT, EFINT,
zspec, zspecFlag, zneur, zneurerr, ztem,z,
j2, av1, err1, zneurFlag, amb2, alb,
alp1, alp2, alp3, alp4,alcirr, alsb, ala220, alagn, alir,
nirtem, errir3,
als12, als25, als60, als90, als100, als110, als140, als160, als250, als350,
als500, als850,
als1250, als1380,
nirflag,
nedtp, sdsstp, nameNED, nameSDSS, name2MASS

----------------------------------------------------------------------------
--------------
The format is

fmt='(a12,1x,2(f10.5,1x),i3,1x,4(f9.3,1x),4(i2,1x),3x,16(f7.2,1x),i3,3x,f13.5,
1x,f10.5,1x,f10.6,1x,i3,1x,4(f10.6,1x),i2,1x,f5.2,1x,f10.3,1x,i3,1x,f7.2,1x,f7.2,3
x,4(f5.2,1x),4(f7.2,1x),3x,f7.2,1x,i4,3x,f10.3,1x,14(f6.2,1x),i3,1x,2(a6,1x),a23,1
x,a22,1x,a22)'
----------------------------------------------------------------------------
-----------------

nameIRAS = IRAS FSC name
ra,dec = positions (priority as for posFlag)
posFlag = position flag (1=SDSS,2=2MASS,3=NVSS,4=NED,5=IRAS FSC, prioritised
in that order: 15=6DF, prioritised after 2MASS, and before NVSS)
s12,s24,s60,s100 are IRAS FSS fluxes (or upper lims) in Jy
nq1, nq2, nq3, nq4 are IRAS flux-quality flags (1=upper limit, 2=low quality, 3=good quality)
am1,am2,am3,am4,am5,am6,am7,am8 are ugriJHK magnitudes
emi, i=1-8, are corresponding errors
photFlag = photometry flag (1=2MASS XSC, 2=2MASS PSC)
FINT, EFINT are 1.4 GHz flux and error
zspec = spectroscopic redshift
zspecFlag = spec z flag (0=unavailable, 1=SDSS, 2=PSCz, 3=FSSz, 4=6dF,
5=NED)
zneur = neural network photometric z
zneurerr = neural network photometric z error
ztem = template-fitting z
z = adopted redshift
j2= optical galaxy template type (1=E, 2=Sab, 3=Sbc, 4=Scd, 5=Sdm, 6=sb, 7=QSO1, 8=QSO2,
    9=QSO3, see Rowan-Robinson et al 2008 for details)
av1= Av from optical galaxy template fit
err1= reduced chi^2 for galaxy template fit
zneurFlag = neural network z flag (1=2MASS, 2=SDSS, 3=NVSS, prioritised in
the same order)
amb2= MB for adopted z
alb= optical bolometric luminosity in solar units
alp1= fraction of contribution at 60 mu of cirrus ir template
alp2= fraction of contribution at 60 mu of M82 ir template
alp3= fraction of contribution at 60 mu of AGN torus ir template
alp4= fraction of contrinution at 60 mu of A220 ir template
alcirr= bolometric luminosity in cirrus component, in solar units
alsb= bolometric luminosity in M82 component
alagn= bolometric luminosity in AGN dust torus component
ala220= bolometric luminosity in A220 component
nirtem= ir template type (1= cirrus dominated, 2= M82 dominated, 3= A220 dominated
     4= AGN dust torus dominated, at 60 mu, 5= detected at 60 mu only)
errir3= reduced chi^2 for ir template fit
als12 ... als1380 are log10(predicted flux in Jy) at 12, 25, 60, 90, 100, 110,140,
    160, 250, 350, 500, 850, 1250, 1380 microns
nirflag = 1 if ir parameters calculated using spectroscopic redshift, =2 if using
neural network redshift, = 3 if using template redshift
nedtp is the NED type flag
sdsstp is the SDSS type flag
nameNED is the NED name
nameSDSS is the SDSS name
name2MASS is the 2MASS name

MRR+LW 18/04/09
