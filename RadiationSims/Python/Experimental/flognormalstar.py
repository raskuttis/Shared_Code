from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_pdf import *
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

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_SDs/'
pdflist = ['sdmasspdfallcirc.dat', 'sdpdfallcirc.dat']
nds = len(pdflist)

plotdens = []
plotpdfs = []
plotpdfms = []
plotfits = []
plotfitms = []

for i in xrange(0,nds):

    pdffile = pdflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    msol = out_msol(outlines)

    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    eff = hst_eff(hstdata, mcloud)

    pdftime = sim_pdftime(pdflines, tff)

    tmini = np.abs(time - 0.1).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf(pdflines,tminpdfi)
    print np.sum(plotpdf)
    exit()
    
    # Calculate the total mass as a consistency check
    # mtot = sim_pdf_mass(pdflines,tminpdfi,(4*rcloud)**2,msol)
    mtot = sim_pdf_mass(pdflines,tminpdfi,4*np.pi*rcloud**2,msol)
    print mtot
    exit()
    
    plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
    plotpdfm = sim_pdfm(pdflines,tminpdfi)
    plotpdfm = np.convolve(plotpdfm, np.ones((nconv,))/nconv, mode='same')
    plotpdf = plotpdf / np.sum(plotpdf)
    nvals = np.sum(plotpdf)
    
    tmintwoi = np.abs(time - 0.9).argmin()
    tmintwo = time[tmintwoi]
    #tmintwo = 0.951
    tminpdftwoi = np.abs(pdftime - tmintwo).argmin()
    
    plotpdftwo = sim_pdf(pdflines,tminpdftwoi)
    plotpdftwo = np.convolve(plotpdftwo, np.ones((nconv,))/nconv, mode='same')
    plotpdftwom = sim_pdfm(pdflines,tminpdftwoi)
    plotpdftwom = np.convolve(plotpdftwom, np.ones((nconv,))/nconv, mode='same')
    plotpdftwo = plotpdftwo / np.sum(plotpdftwo)

    print i, pdffile, tmin, tmintwo
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    plotfit = sim_fitlnhermpdf(plotpdf,plotlogden,plotpdfm,0.15,0.95,nvals)
    plotfitm = sim_fitlnhermpdf(plotpdfm,plotlogden,plotpdfm,0.15,0.95,nvals)
    plotfittwo = sim_fitlnhermpdf(plotpdftwo,plotlogden,plotpdftwom,0.15,0.95,nvals)
    plotfittwom = sim_fitlnhermpdf(plotpdftwom,plotlogden,plotpdftwom,0.15,0.95,nvals)

    plotdens.append(plotden)
    plotpdfs.append(plotpdf)
    plotfits.append(plotfit)
    plotpdfms.append(plotpdfm)
    plotfitms.append(plotfitm)
    plotdens.append(plotden)
    plotpdfs.append(plotpdftwo)
    plotfits.append(plotfittwo)
    plotpdfms.append(plotpdftwom)
    plotfitms.append(plotfittwom)

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
pfid, pall, pmean, pouter = plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[2],plotpdfs[2],'r',plotdens[4],plotpdfs[4],'b',plotdens[6],plotpdfs[6],'g')
plt.plot(plotdens[0],plotfits[0],'--k',plotdens[2],plotfits[2],'--r',plotdens[4],plotfits[4],'--b',plotdens[6],plotfits[6],'--g')
plt.xscale('log')
plt.yscale('log')
plt.legend((pfid, pall), (r"$\displaystyle \Sigma$",r"$\displaystyle \Sigma^{c}_{all}$"))
plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([])

# plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_A$")

plt.subplot(2,2,2)
pfid, pall, pmean, pouter = plt.plot(plotdens[1],plotpdfs[1],'k',plotdens[3],plotpdfs[3],'r',plotdens[5],plotpdfs[5],'b',plotdens[7],plotpdfs[7],'g')
plt.plot(plotdens[1],plotfits[1],'--k',plotdens[3],plotfits[3],'--r',plotdens[5],plotfits[5],'--b',plotdens[7],plotfits[7],'--g')
plt.xscale('log')
plt.yscale('log')
plt.legend((pmean, pouter), (r"$\displaystyle \Sigma^{c}_{*}$",r"$\displaystyle \Sigma^{c}_{max}$"))
plt.text(xcap,ycap,r"$\displaystyle(b)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([])
plt.yticks([])

# plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
# plt.ylabel(r"$\displaystyle P_A$")

plt.subplot(2,2,3)
pfid, pall, pmean, pouter = plt.plot(plotdens[0],plotpdfms[0],'k',plotdens[2],plotpdfms[2],'r',plotdens[4],plotpdfms[4],'b',plotdens[6],plotpdfms[6],'g')
plt.plot(plotdens[0],plotfitms[0],'--k',plotdens[2],plotfitms[2],'--r',plotdens[4],plotfitms[4],'--b',plotdens[6],plotfitms[6],'--g')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(c)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])

plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(2,2,4)
pfid, pall, pmean, pouter = plt.plot(plotdens[1],plotpdfms[1],'k',plotdens[3],plotpdfms[3],'r',plotdens[5],plotpdfms[5],'b',plotdens[7],plotpdfms[7],'g')
plt.plot(plotdens[1],plotfitms[1],'--k',plotdens[3],plotfitms[3],'--r',plotdens[5],plotfitms[5],'--b',plotdens[7],plotfitms[7],'--g')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(d)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])
plt.yticks([])

plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
# plt.ylabel(r"$\displaystyle P_M$")

pp = PdfPages(plotdir + 'flnstar.pdf')
pp.savefig()
pp.close()


