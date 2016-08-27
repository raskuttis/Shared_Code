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

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_3/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_Alt_SDs/'
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
pdftime = sim_pdftime(pdflines, tff)
plotlogden = sim_pdfx(pdflines)
plotden = 10**plotlogden / msol
ndens = np.size(plotden)
npdfts = len(pdftime)
nconv = 25
print 'Finished Reading'

plotdens = []
plotvs = []
plotpdfs = []
nconv = 25
for i in xrange(0,nts):
    
    if ttypelist[i] == 1:
        tmini = np.abs(mstar - tlist[i] * eff).argmin()
    else:
        tmini = np.abs(mof - tlist[i] * (1.0 - eff)).argmin()
    tmin = time[tmini]
    efft = mstar[tmini]
    tmineff = tbreak + (efft - epsbreak) / epsff
    sigmaemin = sigmaedd * (efft) / (efft + 1.0)
    sdintmax = min(3.0, np.log10(sigmaemin))
    sdintmin = -2.0
    sigmaln = 1.5
    sdmean = mcloud / (np.pi * rcloud**2 * msol)
    muln = np.log(sdmean) + 0.5 * sigmaln**2
    scale = np.exp(muln)
    massfrac = lognorm.cdf(10**sdintmax,sigmaln,scale=scale) - lognorm.cdf(10**sdintmin,sigmaln,scale=scale)
    print i, tmin, tmineff, sigmaemin, sdintmax, massfrac
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vmasspdflines,tminpdfi)
    plotpdforig = plotpdf
    plotv = sim_pdfx(vmasspdflines)
    goodplotvs = np.where(plotv > 0.0)
    plotpdf = plotpdf[goodplotvs]
    plotv = plotv[goodplotvs]
    plotvorig = plotv
    mnorm = mcloud / (dx**3)
    plotpdfs.append(np.convolve(plotpdf/mnorm, np.ones((nconv,))/nconv, mode='same'))
    plotvs.append(plotv)



datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_Alt_SDs/'
print 'Reading Out and Hst'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tlist = [1.1,1.88,2.45,3.19]
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
mgas = hst_num(hstdata, 2) / mcloud
mstar = hst_num(hstdata, 13) / mcloud
eff = np.max(mstar)
mof = mgas[0] - mgas - mstar
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
pdftime = sim_pdftime(pdflines, tff)
plotlogden = sim_pdfx(pdflines)
plotden = 10**plotlogden / msol
ndens = np.size(plotden)
npdfts = len(pdftime)
nconv = 25
print 'Finished Reading'

plotdens = []
plotvrs = []
plotpdfrs = []
nconv = 25
for i in xrange(0,nts):
    
    if ttypelist[i] == 1:
        tmini = np.abs(mstar - tlist[i] * eff).argmin()
    else:
        tmini = np.abs(mof - tlist[i] * (1.0 - eff)).argmin()
    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    efft = mstar[tmini]
    tmineff = tbreak + (efft - epsbreak) / epsff
    sigmaemin = sigmaedd * (efft) / (efft + 1.0)
    sdintmax = min(3.0, np.log10(sigmaemin))
    sdintmin = -2.0
    sigmaln = 1.5
    sdmean = mcloud / (np.pi * rcloud**2 * msol)
    muln = np.log(sdmean) + 0.5 * sigmaln**2
    scale = np.exp(muln)
    massfrac = lognorm.cdf(10**sdintmax,sigmaln,scale=scale) - lognorm.cdf(10**sdintmin,sigmaln,scale=scale)
    print i, tmin, tmineff, sigmaemin, sdintmax, massfrac
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vmasspdflines,tminpdfi)
    plotpdforig = plotpdf
    plotv = sim_pdfx(vmasspdflines)
    goodplotvs = np.where(plotv > 0.0)
    plotpdf = plotpdf[goodplotvs]
    plotv = plotv[goodplotvs]
    plotvorig = plotv
    mnorm = mcloud / (dx**3)
    plotpdfrs.append(np.convolve(plotpdf/mnorm, np.ones((nconv,))/nconv, mode='same'))
    plotvrs.append(plotv)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
t1, t2, t3, t4 = plt.plot(plotvs[0],plotpdfs[0],'k',plotvs[1],plotpdfs[1],'r',plotvs[2],plotpdfs[2],'b',plotvs[3],plotpdfs[3],'g')
plt.legend((t1,t2,t3,t4), (r"$\displaystyle t_{\rm 10}$",r"$\displaystyle t_{\rm 50}$",r"$\displaystyle t_{\rm 90}$",r"$\displaystyle t_{\rm of, 50}$"),prop={'size':8})
#plt.xscale('log')
plt.yscale('log')
plt.text(0.05*20,10**(0.9*5-6),r"$\displaystyle(a)$")
plt.axis([0.0,50.0,1.0e-6,1.0e-1])
plt.xticks([0.0,25.0,50.0])
plt.xlabel(r"$\displaystyle \mid v \mid / [{\rm km~s^{-1}}]$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(1,2,2)
plt.plot(plotvrs[0],plotpdfrs[0],'k',plotvrs[1],plotpdfrs[1],'r',plotvrs[2],plotpdfrs[2],'b',plotvrs[3],plotpdfrs[3],'g')
plt.yscale('log')
plt.text(0.05*20,10**(0.9*5-6),r"$\displaystyle(b)$")
plt.axis([0.0,50.0,1.0e-6,1.0e-1])
plt.xticks([0.0,25.0,50.0])
plt.yticks([])
plt.xlabel(r"$\displaystyle \mid v \mid / [{\rm km~s^{-1}}]$")

pp = PdfPages(plotdir + 'ft6-10d.pdf')
pp.savefig()
pp.close()

