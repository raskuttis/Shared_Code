from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_fluxes import *
import numpy as np
import scipy.integrate as spint
import scipy.interpolate as spinterp
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

mlist = ['5.0e3', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '15.0', '20.0', '15.0', '20.0', '25.0', '35.0', '15.0', '25.0', '35.0']
nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflist = np.core.defchararray.add(dflist, ['_N128_Tf4_Fluxes/'] * nfs)

modelsigma = []
modeleff = []
modeleffof = []
tplot = [0.1,0.5,0.9,0.5,3,8]
ttypelist = [1,1,1,0,2,2]
ntplot = len(tplot)

allfabs = np.zeros((ntplot,nfs))
allfabscum = np.zeros((ntplot,nfs))
allsigma = np.zeros(nfs)
allfedd = np.zeros((ntplot,nfs))

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
    
    tff = out_tff(outlines)
    tMyr = out_tMyr(outlines)
    tffMyr = tff / tMyr
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    psi = out_Psi(outlines)
    lsol = psi * msol
    clight = out_clight(outlines)
    kappa = out_kappa(outlines)
    rcloud = out_rcloud(outlines)
    vturb = out_vturb(outlines)
    
    time = hst_time(hstdata, tff)
    eff = hst_mstar(hstdata, mcloud)
    effgas = hst_mgas(hstdata, mcloud)
    effof = effgas[0] - effgas - eff
    eps = hst_eff(hstdata, mcloud)
    tstar = hst_tstar(hstdata, mcloud, tff)
    print i, datafolder, max(time), sigmacloud, tstar
    
    timef, frcl, routs = fluxvst(fluxdata, 6, 2.0, rcloud, tff)
    timef, frclint, routs = fluxintsphervst(fluxdata, 6, 2.0, rcloud, tff)
    fstar = spinterp.interp1d(time, eff * mcloud)
    lstarf = psi * fstar(timef)
    lstarfint = spint.cumtrapz(lstarf, timef * tff, initial=0.0)
    fabsfcl = 1.0 - 4.0 * np.pi * (routs * rcloud)**2 * frcl / lstarf
    fabsfclint = 1.0 - frclint / lstarfint
    lesc = 4.0 * np.pi * (routs * rcloud)**2 * frcl / lsol
    
    allsigma[i] = sigmacloud

    for j in xrange(0,ntplot):
        
        if ttypelist[j] == 2:
            tmini = np.abs(time - tplot[j] / tffMyr - tstar).argmin()
        elif ttypelist[j] == 1:
            tmini = np.abs(eff - tplot[j] * eps).argmin()
        else:
            tmini = np.abs(effof - tplot[j] * (1.0 - eps)).argmin()
        thist = time[tmini]
        tminfi = np.abs(timef - thist).argmin()
        allfabs[j,i] = fabsfcl[tminfi]
        #allfabscum[j,i] = fabsfclint[tminfi]
        allfabscum[j,i] = lesc[tminfi]
    
        rad, rhofrint = fluxintoutvsr(fluxdata, 7, thist, rcloud, tff)
        rad, rhodphidrint = fluxintoutvsr(fluxdata, 15, thist, rcloud, tff)
        feddcum = rhofrint * kappa / (rhodphidrint * clight)
        allfedd[j,i] = feddcum[-1]
        print j, fabsfcl[tminfi], lesc[tminfi]

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
#plt.subplots_adjust(bottom=0.2)
plow, pmid, phigh, pend = plt.plot(allsigma,allfabs[0,:],'k.',allsigma,allfabs[1,:],'r.',allsigma,allfabs[2,:],'b.',allsigma,allfabs[3,:],'g.')
plt.text(10**(1+0.05*2),0.9,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1.0])
plt.yticks([0,0.5,1])
plt.xticks([1.0e1,1.0e2,1.0e3],[' ',' ',' '])

plt.ylabel(r"$\displaystyle f_{\rm abs}$")
#plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

axfedd = plt.subplot(2,2,2)
axfedd.yaxis.set_label_position("right")
plow, pmid, phigh, pend = plt.plot(allsigma,allfabscum[0,:],'k.',allsigma,allfabscum[1,:],'r.',allsigma,allfabscum[2,:],'b.',allsigma,allfabscum[3,:],'g.')
plt.legend((plow,pmid,phigh,pend), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of,final}$"),prop={'size':8},loc=4)
plt.text(10**(1+0.9*2),10**(1+0.9*4),r"$\displaystyle(b)$")
plt.xscale('log')
plt.yscale('log')
plt.axis([1.0e1,1.0e3,1.0e1,1.0e5])
#plt.yticks([])
plt.xticks([1.0e1,1.0e2,1.0e3],[' ',' ',' '])

plt.ylabel(r"$\displaystyle L_{\rm esc} / L_{\odot}$")
#plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(2,2,3)
plt.subplots_adjust(bottom=0.2)
plow, pmid = plt.plot(allsigma,allfabs[4,:],'k.',allsigma,allfabs[5,:],'r.')
plt.text(10**(1+0.05*2),0.9,r"$\displaystyle(c)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1.0])
plt.yticks([0,0.5,1])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle f_{\rm abs}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

axfedda = plt.subplot(2,2,4)
axfedda.yaxis.set_label_position("right")
plow, pmid = plt.plot(allsigma,allfabscum[4,:],'k.',allsigma,allfabscum[5,:],'r.')
plt.legend((plow,pmid), (r"$\displaystyle t = t_* + 3~{\rm Myr}$",r"$\displaystyle t = t_* + 8~{\rm Myr}$"),prop={'size':8},loc=4)
plt.text(10**(1+0.9*2),10**(2.69+0.9*2.31),r"$\displaystyle(d)$")
plt.xscale('log')
plt.yscale('log')
plt.axis([1.0e1,1.0e3,5.0e2,1.0e5])
#plt.yticks([])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle L_{\rm esc} / L_{\odot}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'ffabssigma.pdf')
pp.savefig()
pp.close()

