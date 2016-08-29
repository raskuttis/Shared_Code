from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_fluxes import *
from ..Hyperion.hyp_models import *
import numpy as np
import scipy.integrate as spint
import scipy.interpolate as spinterp
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B50.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B20.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B5.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B2.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B1.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.5_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.2_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.1_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_Fluxes/']
betas = [50.0, 20.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]
nfs = len(dflist)

modelsigma = []
modeleff = []
modeleffof = []
tplot = [0.1,0.5,0.9,0.5,3,8]
ttypelist = [1,1,1,0,2,2]
ntplot = len(tplot)
bzwrong = 1

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
    
    cs = out_csound(outlines)
    gravc = out_G(outlines)
    rhobar = out_rhocloud(outlines)
    beta = betas[i]
    bz = np.sqrt(8.0 * np.pi * cs**2 * rhobar / beta)
    if bzwrong == 1:
        bz = bz * np.sqrt(4.0 * np.pi)
    mf0 = 2 * np.sqrt(gravc) * mcloud / (bz * rcloud**2)
    allsigma[i] = mf0
    
    time = hst_time(hstdata, tff)
    eff = hst_num(hstdata, 13) / mcloud
    effgas = hst_num(hstdata, 2) / mcloud
    effof = effgas[0] - effgas - eff
    eps = np.max(eff)
    tmini = np.min(np.where(eff > 0.0))
    tstar = time[tmini]
    print i, datafolder, max(time), sigmacloud, tstar
    
    #allsigma[i] = sigmacloud

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
plt.legend((plow, pmid, phigh, pend), (r"$\displaystyle t_{\rm 10}$",r"$\displaystyle t_{\rm 50}$",r"$\displaystyle t_{\rm 90}$",r"$\displaystyle t_{\rm of, 50}$"),prop={'size':8},loc=3)
plt.text(10**(-0.3+0.05*2.3),0.9,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([0.5,1.0e2,-0.2,1.0])
plt.yticks([0,0.5,1])
plt.xticks([1.0,1.0e1,1.0e2])

plt.ylabel(r"$\displaystyle f_{\rm abs, F}$")
plt.xlabel(r"$\displaystyle \mu_{\Phi,0}$")

plt.subplot(1,2,2)
plt.subplots_adjust(bottom=0.2)
plow, pmid = plt.plot(allsigma,allfabs[4,:],'k.',allsigma,allfabs[5,:],'r.')
plt.legend((plow,pmid), (r"$\displaystyle t = t_* + 3~{\rm Myr}$",r"$\displaystyle t = t_* + 8~{\rm Myr}$"),prop={'size':8},loc=3)
plt.text(10**(-0.3+0.05*2.3),0.9,r"$\displaystyle(b)$")
plt.xscale('log')
plt.axis([0.5,1.0e2,-0.2,1.0])
plt.yticks([0,0.5,1],[' ',' ',' '])
plt.xticks([1.0,1.0e1,1.0e2])

#plt.ylabel(r"$\displaystyle f_{\rm abs}$")
plt.xlabel(r"$\displaystyle \mu_{\Phi,0}$")

pp = PdfPages(plotdir + 'ft6-10a.pdf')
pp.savefig()
pp.close()

