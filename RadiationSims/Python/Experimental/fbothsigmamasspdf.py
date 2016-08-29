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

tlist = [0.1,0.5,0.9,0.98]
nts = len(tlist)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
clight = out_clight(outlines)
rcloud = out_rcloud(outlines)
msol = out_msol(outlines)
dx = out_dx(outlines)
vturb = out_vturb(outlines)

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

    tmini = np.abs(mstar - tlist[i] * eff).argmin()
    tmin = time[tmini]
    sigmaedd = sigmaeddt[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    plotpdfs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same'))
    
    fitpdf = plotpdf[np.isfinite(plotpdf)]
    fitden = plotlogden[np.isfinite(plotpdf)]
    sigmaxi = np.argmax(fitpdf)
    pouts = np.polyfit(fitden[sigmaxi:], fitpdf[sigmaxi:], 2)
    plotvth = pouts[0] * plotlogden**2 + pouts[1] * plotlogden + pouts[2]
    plotdens.append(plotden)
    plotvths.append(plotvth)

pdffile = 'sdmasspdfallcirc.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'rsdmasspdfallcirc.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdenrs = []
plotpdfrs = []
plotvrths = []
plotrnorms = []
nconv = 25
for i in xrange(0,nts):
    
    tmini = np.abs(mstar - tlist[i] * eff).argmin()
    tmin = time[tmini]
    sigmaedd = sigmaeddt[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    plotpdfrs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same'))
    
    fitpdf = plotpdf[np.isfinite(plotpdf)]
    fitden = plotlogden[np.isfinite(plotpdf)]
    sigmaxi = np.argmax(fitpdf)
    pouts = np.polyfit(fitden[sigmaxi:], fitpdf[sigmaxi:], 2)
    plotvth = pouts[0] * plotlogden**2 + pouts[1] * plotlogden + pouts[2]
    plotdenrs.append(plotden)
    plotvrths.append(plotvth)
    plotrnorms.append(plotpdf**2 / rcloud**2)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
t1, t2, t3, t4 = plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[1],plotpdfs[1],'r',plotdens[2],plotpdfs[2],'b',plotdens[3],plotpdfs[3],'g')
plt.plot(plotdens[0],plotvths[0],'--k',plotdens[1],plotvths[1],'--r',plotdens[2],plotvths[2],'--b',plotdens[3],plotvths[3],'--g')
plt.legend((t1,t2,t3,t4), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.99 \varepsilon_{\rm final}$"),prop={'size':8})
plt.xscale('log')
#plt.yscale('log')
plt.text(10**(0.05*6-2),0.9*25,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,0.0,25.0])
#plt.xticks([])
plt.ylabel(r"$\displaystyle v_r / {\rm km s^{-1}}$")
plt.xlabel(r"$\displaystyle \Sigma^c / M_{\odot}~{\rm pc^{-2}}$")
plt.yticks([0.0,10.0,20.0])

plt.subplot(1,2,2)
plt.plot(plotdens[0] / plotrnorms[0],plotpdfs[0],'k',plotdens[1] / plotrnorms[1],plotpdfs[1],'r',plotdens[2] / plotrnorms[2],plotpdfs[2],'b',plotdens[3] / plotrnorms[3],plotpdfs[3],'g')
plt.xscale('log')
#plt.yscale('log')
plt.text(10**(0.05*6-2),0.9*25,r"$\displaystyle(b)$")
plt.axis([1.0e-2,1.0e4,0.0,25.0])
#plt.xticks([])
plt.yticks([0.0,10.0,20.0],[' ',' ',' '])
plt.xlabel(r"$\displaystyle \Sigma^c / M_{\odot}~{\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'fvsigma.pdf')
pp.savefig()
pp.close()

