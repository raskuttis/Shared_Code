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

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_Alt_SDs/'
print 'Reading Out and Hst'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tlist = [0.1,0.5,0.9,0.5]
ttypelist = [1,1,1,0]
nts = len(tlist)

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
mgas = hst_mgas(hstdata, mcloud)
mof = mgas[0] - mgas - mstar
eff = hst_eff(hstdata, mcloud)
psi = out_Psi(outlines)
psicgs = 2000.0

print 'Reading Densities'
sigmaedd = surf_edd(psicgs)
epsff = 0.3
tbreak = 0.6
epsbreak = 0.0
tbrmax = 1.9
vnorm = gravc * mcloud / (4.0 * rcloud**2) * tff * epsff
pdffile = 'sdmasspdfallcirc.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vsdmasspdfallcirc.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vlrmasspdfall.dat'
vmasspdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vrlrmasspdfall.dat'
vrmasspdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
print 'Reading Log Densities'
pdffile = 'vlogmasspdfall.dat'
vlogmasspdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vrlogmasspdfall.dat'
vrlogmasspdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotlogden = sim_pdfx(pdflines)
plotden = 10**plotlogden / msol
ndens = np.size(plotden)
npdfts = len(pdftime)
nconv = 25
print 'Finished Reading'

allvs = np.zeros((ndens,npdfts))
allrs = np.zeros((ndens,npdfts))
allsds = np.zeros((ndens,npdfts))
for i in xrange(0, ndens):
    sigmain = plotden[i]
    calctime = np.asarray(pdftime)
    rout, vout, tmax = vsigmat_out(sigmain, calctime, psicgs, epsff, tbreak, epsbreak, tbrmax, outlines)
    rnorm = rout / rcloud
    allvs[i,:] = vout
    allrs[i,:] = rnorm
    allsds[i,:] = sigmain / rnorm**2

ndenhr = 10000
nplotts = 10000
plotdenhr = np.logspace(-2.0,3.0,num=ndenhr)
plottimes = np.linspace(0.0,10.0,num=nplotts)
epsstarts = plotdenhr / sigmaedd / (1.0 - plotdenhr / sigmaedd)
allvhrs = np.zeros((ndenhr,npdfts))
allvinf = np.zeros((ndenhr,nplotts))
allrhrs = np.zeros((ndenhr,npdfts))
allsdhrs = np.zeros((ndenhr,npdfts))
for i in xrange(0, ndenhr):
    sigmain = plotdenhr[i]
    calctime = np.asarray(pdftime)
    rout, vout, tmax = vsigmat_out(sigmain, calctime, psicgs, epsff, tbreak, epsbreak, tbrmax, outlines, rmfac=2.0)
    rnorm = rout / rcloud
    allvhrs[i,:] = vout
    allrhrs[i,:] = rnorm
    allsdhrs[i,:] = sigmain / rnorm**2
    rmax = sigmain * msol * kappa
    #print i, sigmain, rmax
    calctime = plottimes
    rout, vout, tmax = vsigmat_out(sigmain, calctime, psicgs, epsff, tbreak, epsbreak, tbrmax, outlines, rmfac=rmax)
    allvinf[i,:] = vout

plotdens = []
plotvs = []
plotlogvs = []
plotpdfs = []
plotpdfrs = []
plotlogpdfs = []
plotlogpdfrs = []
plotpdfths = []
plotpdfrths = []
plotvrs = []
plotlogvrs = []
plotvths = []
plotdenths = []
plotvinfs = []
plotpdfinfs = []
plotpdfrinfs = []
nconv = 25
for i in xrange(0,nts):
    
    if ttypelist[i] == 1:
        tmini = np.abs(mstar - tlist[i] * eff).argmin()
    else:
        tmini = np.abs(mof - tlist[i] * (1.0 - eff)).argmin()
    tmin = time[tmini]
    tmineff = tbreak + (tlist[i] * eff - epsbreak) / epsff
    sigmaemin = sigmaedd * (tlist[i] * eff) / (tlist[i] * eff + 1.0)
    sdintmax = min(3.0, np.log10(sigmaemin))
    sdintmin = -2.0
    sigmaln = 0.85
    sdmean = mcloud / (np.pi * rcloud**2 * msol)
    muln = np.log(sdmean) + 0.5 * sigmaln**2
    scale = np.exp(muln)
    massfrac = lognorm.cdf(10**sdintmax,sigmaln,scale=scale) - lognorm.cdf(10**sdintmin,sigmaln,scale=scale)
    print i, tmin, tmineff, sigmaemin, sdintmax, massfrac
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vlogmasspdflines,tminpdfi)
    plotpdflogorig = plotpdf
    plotlogv = sim_pdfx(vlogmasspdflines)
    plotv = 10**plotlogv
    goodplotvs = np.where(plotv > 0.0)
    plotpdf = plotpdf[goodplotvs]
    plotv = plotv[goodplotvs]
    plotvorig = plotv
    mnorm = mcloud / (dx**3)
    plotlogpdfs.append(np.convolve(plotpdf/mnorm, np.ones((nconv,))/nconv, mode='same'))
    plotlogvs.append(plotv)
    
    plotpdf = sim_pdf_zip(vrlogmasspdflines,tminpdfi)
    plotlogv = sim_pdfx(vrlogmasspdflines)
    plotv = 10**plotlogv
    goodplotvs = np.where(plotv > 0.0)
    plotpdf = plotpdf[goodplotvs]
    plotv = plotv[goodplotvs]
    mnorm = mcloud / (dx**3)
    plotlogpdfrs.append(np.convolve(plotpdf/mnorm, np.ones((nconv,))/nconv, mode='same'))
    plotlogvrs.append(plotv)
    plotlogvrs[i]
    plotlogpdfrs[i]
    
    plotpdf = sim_pdf_zip(vmasspdflines,tminpdfi)
    plotpdforig = plotpdf
    plotv = sim_pdfx(vmasspdflines)
    goodplotvs = np.where(plotv > 0.0)
    plotpdf = plotpdf[goodplotvs] * plotv[goodplotvs]
    plotv = plotv[goodplotvs]
    plotvorig = plotv
    mnorm = mcloud / (dx**3)
    plotpdfs.append(np.convolve(plotpdf/mnorm, np.ones((nconv,))/nconv, mode='same'))
    plotvs.append(plotv)
    
    plotpdf = sim_pdf_zip(vrmasspdflines,tminpdfi)
    plotv = sim_pdfx(vrmasspdflines)
    goodplotvs = np.where(plotv > 0.0)
    plotpdf = plotpdf[goodplotvs] * plotv[goodplotvs]
    plotv = plotv[goodplotvs]
    mnorm = mcloud / (dx**3)
    plotpdfrs.append(np.convolve(plotpdf/mnorm, np.ones((nconv,))/nconv, mode='same'))
    plotvrs.append(plotv)
    
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    pdenfunc = interp1d(plotlogden,plotdenpdf)
    dlnsigma = plotden[1] - plotden[0]
    
    vout = allvs[:,tminpdfi]
    sigmaout = allsds[:,tminpdfi]
    
    nfac = 2.0
    plotpdfth = np.zeros(np.size(plotvorig))
    plotpdfrth = np.zeros(np.size(plotv))
    for j in xrange(0,ndenhr):
        tsd = plotdenhr[j] * nfac**2
        sdind = np.abs(sigmaout - tsd).argmin()
        tv = vout[sdind]
        sdmean = mcloud / (np.pi * rcloud**2 * msol)
        muln = np.log(sdmean) + 0.5 * sigmaln**2
        scale = np.exp(muln)
        pden = lognorm.pdf(tsd, sigmaln, scale = scale)
        pden = pdenfunc(np.log10(tsd * msol))
        sigmap = tsd / sigmaedd
        sigkappa = kappa * tsd * msol
        if (sigmap > eff / (1.0 + eff)):
            epsstart = eff
        else:
            epsstart = sigmap / (1.0 - sigmap)
        epsav = 0.5 * (epsstart + eff)
        voutsq = gravc * mcloud * epsav / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))
        voutfacsq = gravc * mcloud * epsav / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa) / nfac) - nfac * (1.0 - np.exp(-1.0 * sigkappa / nfac**2)))
        tv = np.sqrt(voutsq - voutfacsq)
        vind = np.abs(plotv - tv).argmin()
        if (vind < 0):
            vind = 0
        if (vind >= np.size(plotv)):
            vind = np.size(plotv) - 1
        else:
            plotpdfrth[vind] = plotpdfrth[vind] + pden
        vind = np.abs(plotvorig - tv).argmin()
        if (vind < 0):
            vind = 0
        if (vind >= np.size(plotvorig)):
            vind = np.size(plotvorig) - 1
        else:
            plotpdfth[vind] = plotpdfth[vind] + pden
        #print j, tsd, tv, vind, pden

    vnorm = np.sum(plotpdfth) / (np.sum(plotpdforig/mnorm) * massfrac)
    vrnorm = np.sum(plotpdfrth) / (np.sum(plotpdf/mnorm) * massfrac)

    #plotpdfths.append(np.convolve(plotpdfth/vnorm, np.ones((nconv,))/nconv, mode='same'))
    #plotpdfrths.append(np.convolve(plotpdfrth/vrnorm, np.ones((nconv,))/nconv, mode='same'))

    vout = allvhrs[:,tminpdfi]
    plotpdfth = np.zeros(np.size(plotvorig))
    plotpdfrth = np.zeros(np.size(plotv))
    for j in xrange(0,ndenhr):
        tsd = plotdenhr[j]
        sigmap = tsd / sigmaedd
        sigkappa = tsd * msol * kappa
        sdmean = mcloud * (1.0 - epsstarts[j]) / (np.pi * rcloud**2 * msol)
        pden = lognorm.pdf(tsd, sigmaln, scale = scale)
        tv = vout[j]
        vind = np.abs(plotv - tv).argmin()
        if (tv < min(plotv)):
            vind = 0
        elif (tv > max(plotv)):
            vind = np.size(plotv) - 1
        else:
            plotpdfrth[vind] = plotpdfrth[vind] + pden
        vind = np.abs(plotvorig - tv).argmin()
        if (tv < min(plotvorig)):
            vind = 0
        elif (tv > max(plotvorig)):
            vind = np.size(plotvorig) - 1
        else:
            plotpdfth[vind] = plotpdfth[vind] + pden

    plotpdfth = plotpdfth * plotvorig
    plotpdfrth = plotpdfrth * plotv
    vnorm = np.sum(plotpdfth) / ((1.0 - eff) * massfrac)
    vrnorm = np.sum(plotpdfrth) / ((1.0 - eff) * massfrac)

    plotpdfths.append(np.convolve(plotpdfth/vnorm, np.ones((nconv,))/nconv, mode='same'))
    plotpdfrths.append(np.convolve(plotpdfrth/vrnorm, np.ones((nconv,))/nconv, mode='same'))

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
t1, t2, t3, t4 = plt.plot(plotvs[0],plotpdfs[0],'k',plotvs[1],plotpdfs[1],'r',plotvs[2],plotpdfs[2],'b',plotvs[3],plotpdfs[3],'g')
#plt.plot(plotlogvs[2],plotlogpdfs[2] + 1.0e-10,'--b',plotlogvs[3],plotlogpdfs[3] + 1.0e-10,'--g')
plt.legend((t1,t2,t3,t4), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of, final}$"),prop={'size':8},loc=3)
plt.xscale('log')
plt.yscale('log')
plt.text(10**(0.05*3-1),10**(0.9*5-5),r"$\displaystyle(a)$")
plt.axis([1.0e-1,1.0e2,1.0e-6,1.0])
plt.xlabel(r"$\displaystyle \mid v \mid / {\rm km s^{-1}}$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(1,2,2)
plt.plot(plotvrs[0],plotpdfrs[0],'k',plotvrs[1],plotpdfrs[1],'r',plotvrs[2],plotpdfrs[2],'b',plotvrs[3],plotpdfrs[3],'g')
#plt.plot(plotlogvrs[2],plotlogpdfrs[2] + 1.0e-10,'--b')
plt.xscale('log')
plt.yscale('log')
plt.text(10**(0.05*3-1),10**(0.9*5-5),r"$\displaystyle(b)$")
plt.axis([1.0e-1,1.0e2,1.0e-6,1.0])
plt.yticks([])
plt.xlabel(r"$\displaystyle v_r / {\rm km s^{-1}}$")

pp = PdfPages(plotdir + 'fvloglogfit.pdf')
pp.savefig()
pp.close()

