from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_fluxes import *
from ..Hyperion.hyp_math import *
from ..Hyperion.hyp_models import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.special import erf
from scipy.stats import lognorm

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r0.dat'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = set_model_lookup('Fluxes/')
nfs = len(dflist)

modelsigma = []
modeleff = []
modeleffof = []
psicgs = 2000.0
sigmaedd = surf_edd(psicgs)

altfabs = np.zeros(nfs)
allfabs = np.zeros(nfs)
allsigma = np.zeros(nfs)
allfedd = np.zeros((3,nfs))
allvinf = np.zeros((nfs, 3))
altvinf = np.zeros((nfs, 3))
alleps = np.zeros(nfs)

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    psi = out_Psi(outlines)
    clight = out_clight(outlines)
    kappa = out_kappa(outlines)
    rcloud = out_rcloud(outlines)
    vturb = out_vturb(outlines)
    vesc = vturb * np.sqrt(5.0 * 2.0 / (3.0 * 2.0))
    gravc = out_G(outlines)
    eps = hst_eff(hstdata, mcloud)
    msol = out_msol(outlines)
    
    time = hst_time(hstdata, tff)
    eff = hst_mstar(hstdata, mcloud)
    eps = hst_eff(hstdata, mcloud)
    alleps[i] = eps
    #pnorm = mcloud * vturb
    pnorm = eps * mcloud * vesc
    pnormalt = mcloud * (1.0 - eps) * vesc
    print i, datafolder, max(time), sigmacloud, pnorm
    fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
    
    timefout, prejout = fluxintvst(fluxdata, 12, 2.0, rcloud, tff)
    
    allsigma[i] = sigmacloud
    allfabs[i] = np.max(prejout) / pnorm
    altfabs[i] = np.max(prejout) / pnormalt

## Estimate the final wind velocity based on calculations of v-infinity over the assumed initial distribution of
## surface density and radius
    ndenhr = 1000
    npv = 10000
    plotdenhr = np.logspace(-2.0,4.0,num=ndenhr)
    plotvs = np.linspace(0.0,100.0,num=npv)
    sigmaln = 1.0
    sdmean = mcloud / (np.pi * rcloud**2 * msol)
    muln = np.log(sdmean) + 0.5 * sigmaln**2
    scale = np.exp(muln)
    plotpdfth = np.zeros(npv)
    plotpdfthalt = np.zeros(npv)
    njkrs = 100
    rlist = np.linspace(0.1, 2.0, num = njkrs)
    probr = rlist
    probr = probr / np.sum(probr)
    for j in xrange(0,ndenhr):
        tsd = plotdenhr[j]
        pden = lognorm.pdf(tsd, sigmaln, scale = scale)
        sigmap = tsd / sigmaedd
        sigkappa = kappa * tsd * msol
        epsav = eps
        for jk in xrange(0,njkrs):
            nrmin = rlist[jk]
            nfac = 2.0 / nrmin
            voutsq = gravc * mcloud * epsav / (rcloud * nrmin * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))
            voutfacsq = gravc * mcloud * epsav / (rcloud * nrmin * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa) / nfac) - nfac * (1.0 - np.exp(-1.0 * sigkappa / nfac**2)))
            tv = np.sqrt(voutsq)
            tvalt = np.sqrt(voutsq - voutfacsq)
            vind = np.abs(plotvs - tv).argmin()
            vindalt = np.abs(plotvs - tvalt).argmin()
            if (vind < 0):
                vind = 0
            if (vind >= npv):
                vind = npv - 1
            else:
                plotpdfth[vind] = plotpdfth[vind] + pden * probr[jk]
            
            if (vindalt < 0):
                vindalt = 0
            if (vindalt >= npv):
                vindalt = npv - 1
            else:
                plotpdfthalt[vindalt] = plotpdfthalt[vindalt] + pden * probr[jk]

    sigmap = sigmacloud / (sigmaedd)
    sigkappa = kappa * sigmacloud * msol
    voutsq = gravc * mcloud * epsav / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))

    plotpdfth = plotpdfth / np.sum(plotpdfth)
    plotpdfthcum = np.cumsum(plotpdfth)
    vthmean = np.sum(plotpdfth * plotvs)
    vmind = np.argmin(np.abs(plotpdfthcum - 0.16))
    vpind = np.argmin(np.abs(plotpdfthcum - 0.84))
    vthm = plotvs[vmind]
    vthp = plotvs[vpind]
    allvinf[i,:] = [(vthmean-vthm) / vesc, vthmean / vesc, (vthp-vthmean) / vesc]

    plotpdfthalt = plotpdfthalt / np.sum(plotpdfthalt)
    plotpdfthcum = np.cumsum(plotpdfthalt)
    vthmean = np.sum(plotpdfthalt * plotvs)
    vmind = np.argmin(np.abs(plotpdfthcum - 0.16))
    vpind = np.argmin(np.abs(plotpdfthcum - 0.84))
    vthm = plotvs[vmind]
    vthp = plotvs[vpind]
    altvinf[i,:] = [(vthmean-vthm) / vesc, vthmean / vesc, (vthp-vthmean) / vesc]
    print i, allfabs[i], altfabs[i], allvinf[i,:], altvinf[i,:], voutsq

    epth = 0.5
    mvals = np.logspace(3.0, 6.0, num=100)
    vuni = np.sqrt(2000.0 * epth / 3.0e10) * (mvals * 2.0e33 * 1000.0 / (1.0 - epth))**(1.0/4.0) / 1.0e5


ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

axv = plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
psim = plt.plot(allsigma,altfabs,'k.')
#plow = plt.errorbar(allsigma,np.squeeze(allvinf[:,1]),yerr=[np.squeeze(allvinf[:,0]),np.squeeze(allvinf[:,2])], fmt='r.')
#pmid = plt.errorbar(allsigma,np.squeeze(altvinf[:,1]),yerr=[np.squeeze(altvinf[:,0]),np.squeeze(altvinf[:,2])], fmt='b.')
#plt.legend((psim, plow, pmid), (r"$\displaystyle \langle p_r \rangle / M$",r"$\displaystyle r = 2r_{\rm cl,0}$",r"$\displaystyle r = \infty$"),prop={'size':8},loc=2)
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,3.0])
plt.yticks([0,1,2,3])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle p_{\rm r} / (M_{\rm of} v_{\rm esc})$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / [M_{\odot}~{\rm pc^{-2}}]$")

pp = PdfPages(plotdir + 'f19.pdf')
pp.savefig()
pp.close()

