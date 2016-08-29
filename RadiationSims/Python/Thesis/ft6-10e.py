from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
from ..Hyperion.hyp_models import *
import matplotlib.pyplot as plt
from ..Hyperion.hyp_math import *
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy.stats import lognorm

## Plot showing velocity distributions

## Define the locations where data is located and where to plot from ..Hyperion.hyp_models
hostname, datadir, hstfile, outfile, plotdir = init_dirs('tiger')
fname = 'ft6-10e.pdf'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_Alt_SDs/'
print 'Reading Out and Hst'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

sdlist = [1.0,10.0,100.0]
nsds = len(sdlist)
tlist = [0.5,0.9,0.5]
ttypelist = [1,1,0]
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

sigmaedd = surf_edd(psicgs)
epsff = 0.25
tbreak = 0.6
epsbreak = 0.0
tbrmax = 1.9
vnorm = gravc * mcloud / (4.0 * rcloud**2) * tff * epsff
pdffile = 'sdmasspdfallcirc.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vsdmasspdfallcirc.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

pdftime = sim_pdftime(pdflines, tff)
plotlogden = sim_pdfx(pdflines)
plotden = 10**plotlogden / msol
ndens = np.size(plotden)
npdfts = len(pdftime)
nconv = 25

plotdenths = []
plotpdfths = []
nconv = 25
for i in xrange(0,nts):
    
    if ttypelist[i] == 1:
        tmini = np.abs(mstar - tlist[i] * eff).argmin()
    else:
        tmini = np.abs(mof - tlist[i] * (1.0 - eff)).argmin()
    tmin = time[tmini]
    print i, tmin
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    plotpdfths.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same'))
    
    plotdenths.append(plotden)

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_Alt_SDs/'
print 'Reading Out and Hst'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

sdlist = [1.0,10.0,100.0]
nsds = len(sdlist)
tlist = [0.5,0.9,0.5]
ttypelist = [1,1,0]
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

sigmaedd = surf_edd(psicgs)
epsff = 0.25
tbreak = 0.6
epsbreak = 0.0
tbrmax = 1.9
vnorm = gravc * mcloud / (4.0 * rcloud**2) * tff * epsff
pdffile = 'sdmasspdfallcirc.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vsdmasspdfallcirc.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

pdftime = sim_pdftime(pdflines, tff)
plotlogden = sim_pdfx(pdflines)
plotden = 10**plotlogden / msol
ndens = np.size(plotden)
npdfts = len(pdftime)
nconv = 25

plotdens = []
plotpdfs = []
nconv = 25
for i in xrange(0,nts):
    
    if ttypelist[i] == 1:
        tmini = np.abs(mstar - tlist[i] * eff).argmin()
    else:
        tmini = np.abs(mof - tlist[i] * (1.0 - eff)).argmin()
    tmin = time[tmini]
    print i, tmin
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    plotpdfs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same'))
    
    plotdens.append(plotden)


ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

axr = plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
t1, t2, t3 = plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[1],plotpdfs[1],'r',plotdens[2],plotpdfs[2],'b')
t1, t2, t3 = plt.plot(plotdenths[0],plotpdfths[0],'--k',plotdenths[1],plotpdfths[1],'--r',plotdenths[2],plotpdfths[2],'--b')
plt.legend((t1,t2,t3), (r"$\displaystyle t_{\rm 50}$",r"$\displaystyle t_{\rm 90}$",r"$\displaystyle t_{\rm of, 50}$"),prop={'size':8})
plt.xscale('log')
plt.axis([0.2,1.0e3,0.0,25.0])
plt.ylabel(r"$\displaystyle v_r / [{\rm km~s^{-1}}]$")
plt.xlabel(r"$\displaystyle \Sigma^c / [M_{\odot}~{\rm pc^{-2}}]$")
plt.yticks([0.0,10.0,20.0])

pp = PdfPages(plotdir + fname)
pp.savefig()
pp.close()

