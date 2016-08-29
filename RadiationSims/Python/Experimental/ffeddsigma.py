from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_fluxes import *
import numpy as np
import scipy.integrate as spint
import scipy.interpolate as spinterp
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r0.dat'
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
dfflist = np.core.defchararray.add(dflist, ['_N128_Tf4_Fluxes/'] * nfs)
dfnflist = np.core.defchararray.add(dflist, ['_N128_Tf4_NF_Fluxes/'] * nfs)

modelsigma = []
modeleff = []
modeleffof = []
tplot = [0.1,0.5,0.9,0.5]
ttypelist = [1,1,1,0]
ntplot = len(tplot)

allfabs = np.zeros((ntplot,nfs))
allfabscum = np.zeros((ntplot,nfs))
allsigma = np.zeros(nfs)
allfedd = np.zeros((ntplot,nfs))
nconv = 20

for i in xrange(0,nfs):
    
    datafolder = dfflist[i]
    datafoldernf = dfnflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
    #fluxnfdata = read_fluxfile(hostname,datadir + datafoldernf + fluxfile)
    
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
    
    allsigma[i] = sigmacloud

    tdel = rcloud / vturb / tff
    for j in xrange(0,ntplot):
        
        if ttypelist[j] == 2:
            tmini = np.abs(time - tplot[j] / tffMyr - tstar).argmin()
            tminidel = np.abs(time - tplot[j] / tffMyr - tstar - tdel).argmin()
        elif ttypelist[j] == 1:
            tmini = np.abs(eff - tplot[j] * eps).argmin()
            tminidel = np.abs(eff - tplot[j] * eps - tdel).argmin()
        else:
            tmini = np.abs(effof - tplot[j] * (1.0 - eps)).argmin()
            tminidel = np.abs(effof - tplot[j] * (1.0 - eps) - tdel).argmin()
        thist = time[tmini]
        rad, rhofr = fluxvsr(fluxdata, 7, thist, rcloud, tff)
        rad, rhodphidr = fluxvsr(fluxdata, 15, thist, rcloud, tff)
        fedd = rhofr * kappa / (rhodphidr * clight)
        fedd = np.convolve(fedd, np.ones((nconv,))/nconv, mode='same')
        goodind = np.argmin(np.fabs(rad - 1.0))
        allfabscum[j,i] = fedd[goodind]
        #allfabs[j,i] = 3.0 - pcoeffsnf[0]
        
        print j, allfabscum[j,i]

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
plow, pmid, phigh, pend = plt.plot(allsigma,allfabscum[0,:],'k.',allsigma,allfabscum[1,:],'r.',allsigma,allfabscum[2,:],'b.',allsigma,allfabscum[3,:],'g.')
plt.legend((plow,pmid,phigh,pend), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of,final}$"),prop={'size':8},loc=1)
#plt.text(10**(1+0.05*2),0.9*2,r"$\displaystyle(a)$")
plt.xscale('log')
plt.yscale('log')
plt.axis([1.0e1,1.0e3,1.0e-1,1.0e3])
#plt.yticks([0,1,2])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle f_{\rm edd}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'ffeddsigma.pdf')
pp.savefig()
pp.close()

