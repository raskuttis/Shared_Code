from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt
from ..Hyperion.hyp_math import *


plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_SDs_Alt/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tlist = [1.0,1.5,2.0,2.5]
sdlist = [0.1,1.0,10.0,100.0]
nts = len(tlist)
nsds = len(sdlist)

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
psicgs = 2000.0

sigmaedd = surf_edd(psicgs)
sigmaeddt = sigmaedd * mstar / (1.0 + mstar)

pdffile = 'sdmasspdfallcirc.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vsdmasspdfallcirc.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdens = []
plotpdfs = []
plotvths = []
nconv = 25
for i in xrange(0,nts):

    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    sigmaedd = sigmaeddt[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    plotvth = np.sqrt(sigmaedd / plotden - 1.) * 0.1 * vturb
    plotpdfs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same'))
    plotdens.append(plotden)
    plotvths.append(plotvth)

sdtimes = []
sdvs = []
sdthtimes = []
sdthvs = []
for i in xrange(0,nsds):

    npdfts = len(pdftime)
    sdind = np.abs(plotden - sdlist[i]).argmin()
    sdtime = []
    sdv = []

    for j in xrange(10,npdfts-10):
        
        plotpdf = sim_pdf_zip(vpdflines,j)
        plotdenpdf = sim_pdf(pdflines,j)
        plotpdf = plotpdf / plotdenpdf
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        sdv.append(plotpdf[sdind])
        sdtime.append(pdftime[j])

    sigmain = sdlist[i] * msol
    vnorm = (psi / (4.0 * np.pi * clight * sigmain) - gravc) * mcloud / rcloud**2 * (1.0 - np.exp(-1.0 * sigmain * kappa))
    eps0 = eps_th(sigmain, psi, gravc, clight)
    tmininds = np.where(mstar > eps0)
    minteg = spint.cumtrapz(mstar[tmininds], time[tmininds], initial=0)
    vout = vnorm * minteg

    sdtimes.append(sdtime)
    sdvs.append(sdv)
    sdthtimes.append(time[tmininds])
    sdthvs.append(vout)


pdffile = 'sdmasspdfmeancirc.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vsdmasspdfmeancirc.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdenrs = []
plotpdfrs = []
plotvrths = []
nconv = 25
for i in xrange(0,nts):
    
    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    sigmaedd = sigmaeddt[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    plotvth = np.sqrt(sigmaedd / plotden - 1.) * 0.1 * vturb
    plotpdfrs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same'))
    plotdenrs.append(plotden)
    plotvrths.append(plotvth)

sdrtimes = []
sdrvs = []
sdrthtimes = []
sdrthvs = []
for i in xrange(0,nsds):
    
    npdfts = len(pdftime)
    sdind = np.abs(plotden - sdlist[i]).argmin()
    sdtime = []
    sdv = []
    
    for j in xrange(10,npdfts-10):
        
        plotpdf = sim_pdf_zip(vpdflines,j)
        plotdenpdf = sim_pdf(pdflines,j)
        plotpdf = plotpdf / plotdenpdf
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        sdv.append(plotpdf[sdind])
        sdtime.append(pdftime[j])

    sigmain = sdlist[i] * msol
    vnorm = (psi / (4.0 * np.pi * clight * sigmain) - gravc) * mcloud / rcloud**2 * (1.0 - np.exp(-1.0 * sigmain * kappa))
    eps0 = eps_th(sigmain, psi, gravc, clight)
    tmininds = np.where(mstar > eps0)
    minteg = spint.cumtrapz(mstar[tmininds], time[tmininds], initial=0)
    vout = vnorm * minteg
    
    sdrtimes.append(sdtime)
    sdrvs.append(sdv)
    sdrthtimes.append(time[tmininds])
    sdrthvs.append(vout)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[1],plotpdfs[1],'r',plotdens[2],plotpdfs[2],'b',plotdens[3],plotpdfs[3],'g')
plt.plot(plotdens[0],plotvths[0],'--k',plotdens[1],plotvths[1],'--r',plotdens[2],plotvths[2],'--b',plotdens[3],plotvths[3],'--g')
plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,0.0,25.0])
#plt.xticks([])

plt.subplot(2,2,2)
plt.plot(plotdenrs[0],plotpdfrs[0],'k',plotdenrs[1],plotpdfrs[1],'r',plotdenrs[2],plotpdfrs[2],'b',plotdenrs[3],plotpdfrs[3],'g')
plt.plot(plotdens[0],plotvrths[0],'--k',plotdens[1],plotvrths[1],'--r',plotdens[2],plotvrths[2],'--b',plotdens[3],plotvrths[3],'--g')
plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,0.0,25.0])
#plt.xticks([])

plt.subplot(2,2,3)
plt.plot(sdtimes[0],sdvs[0],'k',sdtimes[1],sdvs[1],'r',sdtimes[2],sdvs[2],'b',sdtimes[3],sdvs[3],'g')
plt.plot(sdthtimes[0],sdthvs[0],'--k',sdthtimes[1],sdthvs[1],'--r',sdthtimes[2],sdthvs[2],'--b',sdthtimes[3],sdthvs[3],'--g')
#plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([0,3,0.0,25.0])
#plt.xticks([])

plt.subplot(2,2,4)
plt.plot(sdrtimes[0],sdrvs[0],'k',sdrtimes[1],sdrvs[1],'r',sdrtimes[2],sdrvs[2],'b',sdrtimes[3],sdrvs[3],'g')
plt.plot(sdrthtimes[0],sdrthvs[0],'--k',sdrthtimes[1],sdrthvs[1],'--r',sdrthtimes[2],sdrthvs[2],'--b',sdrthtimes[3],sdrthvs[3],'--g')
#plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([0,3,0.0,25.0])
#plt.xticks([])

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()

