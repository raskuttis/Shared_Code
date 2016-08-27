from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_pdf import *
import matplotlib.pyplot as plt
from hyp_math import *
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy.stats import lognorm

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_SDs/'
print 'Reading Out and Hst'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
clight = out_clight(outlines)
rcloud = out_rcloud(outlines)
msol = out_msol(outlines)
dx = out_dx(outlines)
vturb = out_vturb(outlines)
gravc = out_G(outlines)
kappa = out_kappa(outlines)

time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
eff = hst_eff(hstdata, mcloud)
psi = out_Psi(outlines)
tsn = out_tMyr(outlines) * 3.0
psicgs = 2000.0

print 'Reading Densities'
sigmaedd = surf_edd(psicgs)
pdffile = 'sdmasspdfallcirc.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
nconv = 25
print 'Finished Reading'

ndenhr = 10000
npv = 10000
plotdenhr = np.logspace(-2.0,4.0,num=ndenhr)
plotvs = np.linspace(0.0,100.0,num=npv)

sigmaln = 1.0
sdmean = mcloud / (np.pi * rcloud**2 * msol)
muln = np.log(sdmean) + 0.5 * sigmaln**2
scale = np.exp(muln)
    
tmini = np.abs(mstar - 0.1 * eff).argmin()
tmin = time[tmini]
tminpdfi = np.abs(pdftime - tmin).argmin()
    
plotdenpdf = sim_pdf(pdflines,tminpdfi)
plotlogden = sim_pdfx(pdflines)
plotden = 10**plotlogden / msol
pdenfunc = interp1d(plotlogden,plotdenpdf)

plotpdfth = np.zeros(npv)
plotpdfthsim = np.zeros(npv)
njkrs = 100
rlist = np.linspace(0.1, 2.0, num = njkrs)
probr = rlist
probr = probr / np.sum(probr)
for j in xrange(0,ndenhr):
    tsd = plotdenhr[j]
    pden = lognorm.pdf(tsd, sigmaln, scale = scale)
    pdensim = pdenfunc(np.log10(tsd * msol))
    sigmap = tsd / sigmaedd
    sigkappa = kappa * tsd * msol
    epsav = eff
    for jk in xrange(0,njkrs):
        nrmin = rlist[jk]
        voutsq = gravc * mcloud * epsav / (rcloud * nrmin * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))
        tvalt = gravc * mcloud * epsav / (rcloud**2 * nrmin**2 * sigmap) * tsn
        tv = np.sqrt(voutsq)
        tvalt = np.min([tv, tvalt])
        vind = np.abs(plotvs - tv).argmin()
        vindalt = np.abs(plotvs - tvalt).argmin()
        if (vind < 0):
            vind = 0
        if (vind >= npv):
            vind = npv - 1
        else:
            plotpdfth[vind] = plotpdfth[vind] + pden * probr[jk]
            plotpdfthsim[vind] = plotpdfthsim[vind] + pdensim * probr[jk]

vnorm = np.sum(plotpdfth) / (1.0 - eff)
plotpdfth = np.convolve(plotpdfth/vnorm, np.ones((nconv,))/nconv, mode='same')
vthmean = np.sum(plotpdfth * plotvs) / np.sum(plotpdfth)
vnorm = np.sum(plotpdfthsim) / (1.0 - eff)
plotpdfthsim = np.convolve(plotpdfthsim/vnorm, np.ones((nconv,))/nconv, mode='same')
vthsimmean = np.sum(plotpdfthsim * plotvs) / np.sum(plotpdfthsim)

sigmap = plotdenhr / sigmaedd
sigkappa = kappa * plotdenhr * msol
epsav = eff
voutsq = gravc * mcloud * epsav / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))
plotvinfs = np.sqrt(voutsq)

nfac = 2.0
sigmap = plotdenhr / sigmaedd
sigkappa = kappa * plotdenhr * msol
epsstart = sigmap / (1.0 - sigmap)
highinds = np.where(sigmap > eff / (1.0 + eff))
epsstart[highinds] = eff
epsav = 0.5 * (epsstart + eff)
voutsq = gravc * mcloud * epsav / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))
voutfacsq = gravc * mcloud * epsav / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa) / nfac) - nfac * (1.0 - np.exp(-1.0 * sigkappa / nfac**2)))
plotvtwo = np.sqrt(voutsq - voutfacsq)

epsav = eff
print epsav
sigmashell = mcloud * 0.46 / (4.0 * np.pi * rcloud**2 * msol)
sigmap = sigmashell / sigmaedd
sigkappa = kappa * sigmashell * msol
voutsigshell = gravc * mcloud * 0.42 / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))
voutsigshell = np.sqrt(voutsigshell) * 1.25
print vthmean, vthsimmean, voutsigshell

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plt.plot(plotdenhr,plotvinfs,'k',plotdenhr,plotvtwo,'r')
plt.xscale('log')
#plt.yscale('log')
plt.text(10**(0.05*6-2),0.9*50,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,0.0,50.0])
#plt.xticks([])
plt.ylabel(r"$\displaystyle v_r / {\rm km s^{-1}}$")
plt.xlabel(r"$\displaystyle \Sigma^c / M_{\odot}~{\rm pc^{-2}}$")
#plt.yticks([0.0,10.0,20.0],[' ',' ',' '])

plt.subplot(1,2,2)
plt.plot(plotvs,plotpdfth,'k',plotvs,plotpdfthsim,'r')
plt.axvline(vthmean,color='k',linestyle='--')
plt.axvline(vthsimmean,color='r',linestyle='--')
plt.axvline(voutsigshell,color='m',linestyle='--')
#plt.xscale('log')
plt.yscale('log')
plt.text(0.05*50,10**(0.9*3-5),r"$\displaystyle(b)$")
plt.axis([0.0,50.0,1.0e-5,1.0e-2])
plt.yticks([])
plt.xlabel(r"$\displaystyle v_r / {\rm km s^{-1}}$")

pp = PdfPages(plotdir + 'fvinf.pdf')
pp.savefig()
pp.close()

