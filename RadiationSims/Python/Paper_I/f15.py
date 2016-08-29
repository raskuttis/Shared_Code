from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'sdpdfxy.dat'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
nconv = 25

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
msol = out_msol(outlines)

time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
eff = hst_eff(hstdata, mcloud)

pdftime = sim_pdftime(pdflines, tff)

tmini = np.abs(mstar - 0.1 * eff).argmin()
tmin = time[tmini]
tminpdfi = np.abs(pdftime - tmin).argmin()

plotpdf = sim_pdf(pdflines,tminpdfi)
plotpdfm = sim_pdfm(pdflines,tminpdfi)
plotpdf = plotpdf / np.sum(plotpdf)

tmintwoi = np.abs(mstar - 0.5 * eff).argmin()
tmintwo = time[tmintwoi]
tminpdftwoi = np.abs(pdftime - tmintwo).argmin()
plotpdftwo = sim_pdf(pdflines,tminpdftwoi)
plotpdftwom = sim_pdfm(pdflines,tminpdftwoi)
plotpdftwo = plotpdftwo / np.sum(plotpdftwo)

plotlogden = sim_pdfx(pdflines)
plotden = 10**plotlogden / msol
plotfit = sim_fitlnpdf(plotpdf,plotlogden,plotpdfm,0.1,0.9)
plotfitm = sim_fitlnpdf(plotpdfm,plotlogden,plotpdfm,0.1,0.9)
plotfittwo = sim_fitlnpdf(plotpdftwo,plotlogden,plotpdftwom,0.1,0.9)
plotfittwom = sim_fitlnpdf(plotpdftwom,plotlogden,plotpdftwom,0.1,0.9)

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_NF/'
altpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

altpdftime = sim_pdftime(altpdflines, tff)

altplotpdf = sim_pdf(altpdflines,tminpdfi)
altplotpdfm = sim_pdfm(altpdflines,tminpdfi)
altplotpdf = altplotpdf / np.sum(altplotpdf)
altplotpdftwo = sim_pdf(altpdflines,tminpdftwoi)
altplotpdftwom = sim_pdfm(altpdflines,tminpdftwoi)
altplotpdftwo = altplotpdftwo / np.sum(altplotpdftwo)

altplotlogden = sim_pdfx(altpdflines)
altplotden = 10**altplotlogden / msol

altplotfit = sim_fitlnpdf(altplotpdf,plotlogden,altplotpdfm,0.1,0.9)
altplotfitm = sim_fitlnpdf(altplotpdfm,plotlogden,altplotpdfm,0.1,0.9)
altplotfittwo = sim_fitlnpdf(altplotpdftwo,plotlogden,altplotpdftwom,0.1,0.9)
altplotfittwom = sim_fitlnpdf(altplotpdftwom,plotlogden,altplotpdftwom,0.1,0.9)

plotpdfm = np.convolve(plotpdfm, np.ones((nconv,))/nconv, mode='same')
plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
plotpdftwom = np.convolve(plotpdftwom, np.ones((nconv,))/nconv, mode='same')
plotpdftwo = np.convolve(plotpdftwo, np.ones((nconv,))/nconv, mode='same')
altplotpdfm = np.convolve(altplotpdfm, np.ones((nconv,))/nconv, mode='same')
altplotpdf = np.convolve(altplotpdf, np.ones((nconv,))/nconv, mode='same')
altplotpdftwom = np.convolve(altplotpdftwom, np.ones((nconv,))/nconv, mode='same')
altplotpdftwo = np.convolve(altplotpdftwo, np.ones((nconv,))/nconv, mode='same')

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
pfid, pnf = plt.plot(plotden,plotpdf,'k',altplotden,altplotpdf,'r')
plt.plot(plotden,plotfit,'--k',plotden,altplotfit,'--r')
plt.xscale('log')
plt.yscale('log')
plt.legend((pfid, pnf), (r"$\displaystyle {\rm Fiducial}$",r"$\displaystyle {\rm No Feedback}$"))
plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([])

# plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_A$")

plt.subplot(2,2,2)
pfid, pnf = plt.plot(plotden,plotpdftwo,'k',altplotden,altplotpdftwo,'r')
plt.plot(plotden,plotfittwo,'--k',plotden,altplotfittwo,'--r')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(b)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([])
plt.yticks([])

# plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
# plt.ylabel(r"$\displaystyle P_A$")

plt.subplot(2,2,3)
pfid, pnf = plt.plot(plotden,plotpdfm,'k',altplotden,altplotpdfm,'r')
plt.plot(plotden,plotfitm,'--k',plotden,altplotfitm,'--r')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(c)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])

plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(2,2,4)
pfid, pnf = plt.plot(plotden,plotpdftwom,'k',altplotden,altplotpdftwom,'r')
plt.plot(plotden,plotfittwom,'--k',plotden,altplotfittwom,'--r')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(d)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])
plt.yticks([])

plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
# plt.ylabel(r"$\displaystyle P_M$")

pp = PdfPages(plotdir + 'f15.pdf')
pp.savefig()
pp.close()


