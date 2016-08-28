from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_fluxes import *
from hyp_models import *
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

dflist = set_model_lookup('Fluxes/')
nfs = len(dflist)

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
    
    allsigma[i] = sigmacloud

    for j in xrange(0,ntplot):
        
        if ttypelist[j] == 2:
            tmini = np.abs(time - tplot[j] / tffMyr - tstar).argmin()
        elif ttypelist[j] == 1:
            tmini = np.abs(eff - tplot[j] * eps).argmin()
        else:
            tmini = np.abs(effof - tplot[j] * (1.0 - eps)).argmin()
        thist = time[tmini]
        fout = flux_tableouts(fluxdata, thist, 'fabs', rcloud, tff, mcloud, psi, time = time, eff = eff)
        allfabs[j,i] = fout
        
        print j, fout

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
#plt.subplots_adjust(bottom=0.2)
plow, pmid, phigh, pend = plt.plot(allsigma,allfabs[0,:],'k.',allsigma,allfabs[1,:],'r.',allsigma,allfabs[2,:],'b.',allsigma,allfabs[3,:],'g.')
plt.text(10**(1+0.05*2),0.9,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1.0])
plt.yticks([0,0.5,1])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle f_{\rm abs, F}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / [M_{\odot} {\rm pc^{-2}}]$")

plt.subplot(1,2,2)
plt.subplots_adjust(bottom=0.2)
plow, pmid = plt.plot(allsigma,allfabs[4,:],'k.',allsigma,allfabs[5,:],'r.')
plt.text(10**(1+0.05*2),0.9,r"$\displaystyle(b)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1.0])
plt.yticks([0,0.5,1],[' ',' ',' '])
plt.xticks([1.0e1,1.0e2,1.0e3])

#plt.ylabel(r"$\displaystyle f_{\rm abs}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / [M_{\odot} {\rm pc^{-2}}]$")

pp = PdfPages(plotdir + 'f17.pdf')
pp.savefig()
pp.close()

