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

dfflist = set_model_lookup('Fluxes/')
dfnflist = set_model_lookup('NF_Fluxes/')
nfs = len(dfflist)

modelsigma = []
modeleff = []
modeleffof = []
tplot = [0.1,0.5,0.9]
ttypelist = [1,1,1]
ntplot = len(tplot)

allfabs = np.zeros((ntplot,nfs))
allfabscum = np.zeros((ntplot,nfs))
allsigma = np.zeros(nfs)
allfedd = np.zeros((ntplot,nfs))

for i in xrange(2,nfs):
    
    datafolder = dfflist[i]
    datafoldernf = dfnflist[i]
    #print i, datafolder
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
        alpha = flux_tableouts(fluxdata, thist, 'alpha', rcloud, tff, mcloud, 1.0)
        allfabscum[j,i] = alpha
        
        print j, allfabscum[j,i]

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plow, pmid, phigh = plt.plot(allsigma,allfabscum[0,:],'k.',allsigma,allfabscum[1,:],'r.',allsigma,allfabscum[2,:],'b.')
#plt.legend((plow,pmid,phigh), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$"),prop={'size':8},loc=1)
plt.legend((plow,pmid,phigh), (r"$\displaystyle t_{\rm 10}$",r"$\displaystyle t_{\rm 50}$",r"$\displaystyle t_{\rm 90}$"),prop={'size':8},loc=2)
#plt.text(10**(1+0.05*2),0.9*2,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,2.0])
plt.yticks([0,1,2])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle \alpha$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'f5.pdf')
pp.savefig()
pp.close()

