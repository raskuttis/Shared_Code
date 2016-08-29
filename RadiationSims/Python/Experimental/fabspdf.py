from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
nconv = 25

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_SDs/'
pdflist = ['sdpdfallcirc.dat', 'sdpdfallcirc.dat', 'sdpdfmeancirc.dat', 'sdpdfoutercirc.dat']
nds = len(pdflist)

plotdens = []
plotpdfs = []
plotfits = []

for i in xrange(0,nds):

    pdffile = pdflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    msol = out_msol(outlines)
    kappa = out_kappa(outlines)

    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    eff = hst_eff(hstdata, mcloud)

    pdftime = sim_pdftime(pdflines, tff)
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    plottau = 1.0 - np.exp(-1.0 * plotden * kappa)

    tmini = np.abs(mstar - 0.1 * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
    #plotpdf = plotpdf / plotden
    plotpdf = plotpdf / np.sum(plotpdf)
    
    tmintwoi = np.abs(mstar - 0.5 * eff).argmin()
    tmintwo = time[tmintwoi]
    tminpdftwoi = np.abs(pdftime - tmintwo).argmin()

    plotpdftwo = sim_pdf(pdflines,tminpdftwoi)
    plotpdftwo = np.convolve(plotpdftwo, np.ones((nconv,))/nconv, mode='same')
    #plotpdftwo = plotpdftwo / plotden
    plotpdftwo = plotpdftwo / np.sum(plotpdftwo)

    plotdens.append(plottau)
    plotpdfs.append(plotpdf)
    plotdens.append(plottau)
    plotpdfs.append(plotpdftwo)

xcap = 10**(-4.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
pfid, pall, pmean, pouter = plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[2],plotpdfs[2],'r',plotdens[4],plotpdfs[4],'b',plotdens[6],plotpdfs[6],'g')
plt.legend((pfid, pall), (r"$\displaystyle \Sigma$",r"$\displaystyle \Sigma^{c}_{all}$"),prop={'size':8})
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-4,1.0,1.0e-4,1.0e-1])
plt.xticks([1.0e-4,1.0e-2,1.0])

plt.xlabel(r"$\displaystyle f_{\rm abs}$")
plt.ylabel(r"$\displaystyle P_{\Omega}$")

plt.subplot(1,2,2)
pfid, pall, pmean, pouter = plt.plot(plotdens[1],plotpdfs[1],'k',plotdens[3],plotpdfs[3],'r',plotdens[5],plotpdfs[5],'b',plotdens[7],plotpdfs[7],'g')
plt.legend((pmean, pouter), (r"$\displaystyle \Sigma^{c}_{*}$",r"$\displaystyle \Sigma^{c}_{max}$"),prop={'size':8})
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(b)$")
plt.axis([1.0e-4,1.0,1.0e-4,1.0e-1])
plt.xticks([1.0e-4,1.0e-2,1.0])
plt.yticks([])

plt.xlabel(r"$\displaystyle f_{\rm abs}$")
# plt.ylabel(r"$\displaystyle P_M$")

pp = PdfPages(plotdir + 'ffabspdf.pdf')
pp.savefig()
pp.close()


