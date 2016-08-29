from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_fluxes import *
from ..Hyperion.hyp_math import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt
import scipy.interpolate as spinterpol

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hostname = 'raskutti@tiger.princeton.edu'
starfile = 'star'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_NF_Fluxes/'
nffluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_Fluxes/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

msol = out_msol(outlines)
mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
tff = out_tff(outlines)
psi = out_Psi(outlines)
kappa = out_kappa(outlines)

tplot = 1.57

time = hst_time(hstdata, tff)
eff = hst_mstar(hstdata, mcloud)
effgas = hst_mgas(hstdata, mcloud)
stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)

starrad, starmvsr = star_massvsratt(stardata, time, tplot)
starrad = starrad / rcloud
starmvstrcl = star_massvstatr(stardata, 1.0 * rcloud)
starmvstrbox = star_massvstatr(stardata, 2.0 * rcloud)

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_SDs/'
pdffile = 'sdmasspdfallcirc.dat'

pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
nts = len(pdftime)
    
plottime = []
plotfabsth = []
nconv = 25
nconvall = 100

startflag = False
    
for j in xrange(1250,5100):
    
    plotpdf = sim_pdf(pdflines,j)
    plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
    mtot = np.sum(plotpdf) / 1.0e9
    plotpdf = plotpdf / np.sum(plotpdf)
    nvals = np.sum(plotpdf)
    tmini = np.abs(time - pdftime[j]).argmin()
    sigmamean = effgas[tmini] * mcloud / (4.0 * np.pi * rcloud**2)
        
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    if startflag:
        amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
    else:
        amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals)
        
    if (sigma > 0.0):
            
        [ampold, sigmaold, muold] = [amp, sigma, mu]
        sigma = sigma / np.log10(np.exp(1))
        mu = mu / np.log10(np.exp(1))
        mu = mu - sigma**2
        mean = np.exp(mu - 0.5 * sigma**2)
            
        plottime.append(pdftime[j])
        fabsth = fabs_th(sigma,mu,kappa)
        fabsth = 1.0 - np.exp(-1.0 * kappa * sigmamean)
        plotfabsth.append(fabsth)
        startflag = True

plottime = np.asarray(plottime)
plotfabsth = np.asarray(plotfabsth)
#plotfabsth = np.convolve(plotfabsth, np.ones((nconvall,))/nconvall, mode='same')
tminploti = np.abs(plottime - tplot).argmin()
tmini = np.abs(time - tplot).argmin()
mstar = eff[tmini] * mcloud
fabsth = plotfabsth[tminploti]
lstar = psi * mstar

rad, tau = fluxvsr(fluxdata, 3, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
nrads = np.size(rad)
lstarr = np.zeros(nrads)
for i in xrange(0,nrads):
    trad = rad[i]
    starinds = np.where(starrad > trad)
    if (np.size(starinds) > 0):
        starind = np.min(starinds[0])
        lstarr[i] = psi * starmvsr[starind]
    else:
        lstarr[i] = lstar

fabsf = 1.0 - 4.0 * np.pi * (rad * rcloud)**2 * fr / lstarr
fabstau = 1.0 - tau
fabsth = fabsth * rad / rad

timef, taucl, routs = fluxvst(fluxdata, 3, 1.0, rcloud, tff)
fabstaucl = 1.0 - taucl
timef, tauout, routs = fluxvst(fluxdata, 3, 2.0, rcloud, tff)
fabstauout = 1.0 - tauout
fstar = interp1d(time, eff * mcloud)
lstarf = psi * fstar(timef)

timef, frcl, rout = fluxvst(fluxdata, 6, 1.0, rcloud, tff)
fabsfcl = 1.0 - 4.0 * np.pi * (rout * rcloud)**2 * frcl / (psi * starmvstrcl)
timef, frout, rout = fluxvst(fluxdata, 6, 2.0, rcloud, tff)
fabsfout = 1.0 - 4.0 * np.pi * (rout * rcloud)**2 * frout / (psi * starmvstrbox)

timenf, nftauout, nfrouts = fluxvst(nffluxdata, 3, 2.0, rcloud, tff)
fabstaunfout = 1.0 - nftauout

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,ysize])
#plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,2,1)
ft, ff = plt.plot(rad,fabstau,'k', rad,fabsf,'r')
plt.axis([0,2,0,1.0])
plt.xticks([0,1,2])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle f_{\rm abs}$")
plt.text(0.05*2,0.9*1,r"$\displaystyle(a)$")

plt.subplot(1,2,2)
ft, ftnf, ff, fth = plt.plot(timef,fabstauout,'k',timenf,fabstaunfout,'--k',timef,fabsfout,'r',plottime,plotfabsth,'b')
plt.legend((ft, ff, fth), (r"$\displaystyle f_{\rm abs, \tau}$",r"$\displaystyle f_{\rm abs, F}$",r"$\displaystyle f_{\rm abs, th}$"),prop={'size':8},loc=3)
plt.axis([0,2,0,1.0])
plt.xticks([0,1,2])
plt.yticks([0,0.5,1],[' ',' ',' '])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle f_{\rm abs}$")
plt.text(0.05*2,0.9*1,r"$\displaystyle(b)$")

pp = PdfPages(plotdir + 'fflux.pdf')
pp.savefig()
pp.close()




