from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_fluxes import *
from ..Hyperion.hyp_math import *
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

gravc = 6.67e-8
mcloud = 5.0e4 * 2.0e33
rcloud = 15.0 * 3.086e18
psicgs = 2000.0
sigmaedd = surf_edd(psicgs)
ndenhr = 100
npv = 100
epsav = 0.5
fac = gravc * mcloud * epsav / rcloud
kappa = 1000.0
plotdenhr = np.logspace(-2.0,4.0,num=ndenhr)
sigmap = plotdenhr / (sigmaedd)
sigkappa = kappa * plotdenhr * 0.000208908219
voutsq = gravc * mcloud * epsav / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))
vout = np.sqrt(voutsq) / 1.0e5
    
plotvs = np.linspace(0.0,100.0,num=npv)
sigmaln = 2.0
niter = 1000
vouts = np.zeros(ndenhr)
for i in xrange(ndenhr-1,ndenhr):
    muln = np.log(plotdenhr[i]) + 0.5 * sigmaln**2
    logc = np.log10(plotdenhr[i])
    scale = np.exp(muln)
    plotpdfth = np.zeros(npv)
    iterdenhr = np.logspace(logc-4.0,logc+4.0,num=niter)
    medn = lognorm.median(sigmaln, scale=scale)
    pdenm = lognorm.pdf(medn, sigmaln, scale=scale)
    pden = lognorm.pdf(plotdenhr[i], sigmaln, scale = scale)
    print medn, plotdenhr[i], scale, pdenm, pden
    exit()
    for j in xrange(0,niter):
        tsd = iterdenhr[j]
        pden = lognorm.pdf(tsd, sigmaln, scale = scale)
        sigmap = tsd / (sigmaedd)
        sigkappa = kappa * tsd * 0.000208908219
        voutsq = fac / sigmap * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))
        tv = np.sqrt(voutsq) / 1.0e5
        print tsd, pden, tv
        vind = np.abs(plotvs - tv).argmin()
        if (vind < 0):
            vind = 0
        if (vind >= npv):
            vind = npv - 1
        else:
            plotpdfth[vind] = plotpdfth[vind] + pden

    print lognorm.median(sigmaln, scale=scale)
    plotpdfth = plotpdfth / np.sum(plotpdfth)
    vthmean = np.sum(plotpdfth * plotvs)
    vouts[i] = vthmean
    print i, vouts[i]

print np.min(plotpdfth), np.max(plotpdfth)


ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
psim = plt.plot(plotpdfth,plotvs,'k')
plt.xscale('log')
plt.axis([1.0e-6,1.0e4,0,50.0])
#plt.yticks([0,25,50,75,100])
plt.xticks([1.0e-2,1.0e0,1.0e2,1.0e4])

plt.ylabel(r"$\displaystyle v_{\infty} ({\rm km/s})$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot}~{\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'fveffect.pdf')
pp.savefig()
pp.close()

