from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt
from ..Hyperion.hyp_math import *
from scipy.special import erf


plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_SDs/'
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

allvs = np.zeros((ndens,npdfts))
allrs = np.zeros((ndens,npdfts))
allsds = np.zeros((ndens,npdfts))
alltmaxs = np.zeros((ndens,npdfts))
for i in xrange(0, ndens):
    sigmain = plotden[i]
    if (sigmain < sigmaedd):
        calctime = np.asarray(pdftime)
        rout, vout, tmax = vsigmat_out(sigmain, calctime, psicgs, epsff, tbreak, epsbreak, tbrmax, outlines)
        rnorm = rout / rcloud
        allvs[i,:] = vout
        allrs[i,:] = rnorm
        allsds[i,:] = sigmain / rnorm**2
        alltmaxs[i,:] = tmax
    else:
        allvs[i,:] = 0.0
        allrs[i,:] = 1.0
        allsds[i,:] = sigmain
        alltmaxs[i,:] = 1.0e10

sdtimes = []
sdvs = []
sdthtimes = []
sdthvs = []
for i in xrange(0,nsds):
    
    sdind = np.abs(plotden - sdlist[i]).argmin()
    sdtime = []
    sdv = []
    sdthv = []
    
    for j in xrange(1500,npdfts-10):
        
        plotpdf = sim_pdf_zip(vpdflines,j)
        plotdenpdf = sim_pdf(pdflines,j)
        plotpdf = plotpdf / plotdenpdf
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        sdv.append(plotpdf[sdind])
        sdtime.append(pdftime[j])
        thsds = allsds[:,j]
        sdthind = np.abs(thsds - sdlist[i]).argmin()
        #print j, sdthind, thsds[sdthind], sdlist[i], allvs[sdthind,j], allrs[sdthind,j]
        sdthv.append(allvs[sdthind,j])
    
    sdtimes.append(sdtime)
    sdvs.append(sdv)
    sdthtimes.append(sdtime)
    sdthvs.append(sdthv)

plotdens = []
plotpdfs = []
plotvths = []
plotdenths = []
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
    
    nfac = 4.0
    sigmap = plotden / sigmaedd
    sigfrac = sigmap / (epsff * (1.0 - sigmap))
    sigkappa = kappa * plotden * msol
    sigmapfac = sigmap * nfac
    sigkappafac = sigkappa * nfac
    breakfrac = epsbreak / epsff
    calctime = tmin
    
    vout = allvs[:,tminpdfi]
    sigmaout = allsds[:,tminpdfi]
    tmax = alltmaxs[:,tminpdfi]
    epsmax = epsbreak + epsff * (tmax - tbreak)
    highinds = np.where(epsmax > eff)
    epsmax[highinds] = eff
    
    eps0 = sigmapfac / (1.0 - sigmapfac)
    epsav = 0.5 * eps0 + 0.5 * epsmax
    voutsq = gravc * mcloud * epsav / (rcloud * sigmapfac) * (np.sqrt(sigkappafac * np.pi) * erf(np.sqrt(sigkappafac)) - (1.0 - np.exp(-1.0 * sigkappafac)))
    voutfacsq = gravc * mcloud * epsav / (rcloud * sigmapfac) * (np.sqrt(sigkappafac * np.pi) * erf(np.sqrt(nfac * sigkappafac)) - np.sqrt(nfac) * (1.0 - np.exp(-nfac * sigkappafac)))
    #vout = np.sqrt(voutsq - voutfacsq)
    #sigmaout = plotden / nfac
    
    plotdens.append(plotden)
    plotdenths.append(sigmaout)
    plotvths.append(vout)

nfac = 2.0
sigmap = plotden / sigmaedd
sigkappa = kappa * plotden * msol
epsstart = sigmap / (1.0 - sigmap)
highinds = np.where(sigmap > eff / (1.0 + eff))
epsstart[highinds] = eff
epsav = 0.5 * (epsstart + eff)
voutsq = gravc * mcloud * epsav / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa)) - (1.0 - np.exp(-1.0 * sigkappa)))
voutfacsq = gravc * mcloud * epsav / (rcloud * sigmap) * (np.sqrt(sigkappa * np.pi) * erf(np.sqrt(sigkappa) / nfac) - nfac * (1.0 - np.exp(-1.0 * sigkappa / nfac**2)))
plotvinfs = np.sqrt(voutsq - voutfacsq)
plotdinfs = plotden / nfac**2
plotvrinfs = np.zeros(np.size(plotden))
for i in xrange(0, np.size(plotden)):
    sigmain = plotden[i]
    epsstart = sigmain / sigmaedd / (1.0 - sigmain / sigmaedd)
    sigmain = sigmain * (1.0 - epsstart)
    gfrac = 1.0
    divfrac = np.log(eff * sigmaedd / (eff * sigmaedd - 2.0 * gfrac * sigmain))
    rmax = sigmain * msol * kappa / divfrac
    #rmax = 4.0
    calctime = 10.0
    rout, vout, tmax = vsigmat_out(sigmain, calctime, psicgs, epsff, tbreak, epsbreak, tbrmax, outlines, rmfac=rmax, ctmax=1)
    if (plotden[i] > eff / (1.0 + eff) * sigmaedd):
        vout = 0.0
    #print i, sigmain, plotden[i], epsstart, divfrac, rmax, rout, vout
    plotvrinfs[i] = vout

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

axt = plt.subplot(1,2,2)
axt.yaxis.set_label_position("right")
plt.subplots_adjust(bottom=0.2)
t1, t2, t3 = plt.plot(sdtimes[0],sdvs[0],'k',sdtimes[1],sdvs[1],'r',sdtimes[2],sdvs[2],'b')
plt.plot(sdthtimes[0],sdthvs[0],'--k',sdthtimes[1],sdthvs[1],'--r',sdthtimes[2],sdthvs[2],'--b')
plt.legend((t1,t2,t3), (r"$\displaystyle \Sigma^c = 1~M_{\odot}~{\rm pc^{-2}}$",r"$\displaystyle \Sigma^c = 10~M_{\odot}~{\rm pc^{-2}}$",r"$\displaystyle \Sigma^c = 100~M_{\odot}~{\rm pc^{-2}}$"),prop={'size':8},loc=2)
#plt.xscale('log')
#plt.yscale('log')
plt.text(0.9*3,0.9*25,r"$\displaystyle(b)$")
plt.axis([0,3,0.0,25.0])
plt.xticks([0,1,2,3])
plt.ylabel(r"$\displaystyle v_r / {\rm km s^{-1}}$")
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.yticks([0.0,10.0,20.0],[' ',' ',' '])

axr = plt.subplot(1,2,1)
t1, t2, t3 = plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[1],plotpdfs[1],'r',plotdens[2],plotpdfs[2],'b')
plt.plot(plotdenths[0],plotvths[0],'--k',plotdenths[1],plotvths[1],'--r',plotdenths[2],plotvths[2],'--b',plotdinfs,plotvinfs,'m')
plt.legend((t1,t2,t3), (r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of, final}$"),prop={'size':8})
plt.xscale('log')
#plt.yscale('log')
plt.text(10**(0.05*6-2),0.9*25,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,0.0,25.0])
#plt.xticks([])
plt.ylabel(r"$\displaystyle v_r / {\rm km s^{-1}}$")
plt.xlabel(r"$\displaystyle \Sigma^c / M_{\odot}~{\rm pc^{-2}}$")
plt.yticks([0.0,10.0,20.0])

pp = PdfPages(plotdir + 'fvtsigma.pdf')
pp.savefig()
pp.close()

