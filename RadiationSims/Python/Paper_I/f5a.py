from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'denpdfall.dat'
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

tmini = np.argmax(mstar > 0.1 * eff)
tmin = time[tmini]
tminpdfi = np.abs(pdftime - tmin).argmin()
print tmini, tmin, tminpdfi

plotpdf = sim_pdf(pdflines,tminpdfi)
plotpdfm = sim_pdfm(pdflines,tminpdfi)
plotpdfm = np.convolve(plotpdfm, np.ones((nconv,))/nconv, mode='same')
plotpdf = plotpdf / np.sum(plotpdf)

tmintwoi = np.argmax(mstar > 0.5 * eff)
tmintwo = time[tmintwoi]
tminpdftwoi = np.abs(pdftime - tmintwo).argmin()
print tmintwoi, tmintwo, tminpdftwoi
plotpdftwo = sim_pdf(pdflines,tminpdftwoi)
plotpdftwom = sim_pdfm(pdflines,tminpdftwoi)
plotpdftwom = np.convolve(plotpdftwom, np.ones((nconv,))/nconv, mode='same')
plotpdftwo = plotpdftwo / np.sum(plotpdftwo)

plotlogden = sim_pdfx(pdflines)
plotden = 10**plotlogden
plotfit = sim_fitlnpdf(plotpdf,plotlogden,plotpdfm,0.1,0.9)
plotfitm = sim_fitlnpdf(plotpdfm,plotlogden,plotpdfm,0.1,0.9)
plotfittwo = sim_fitlnpdf(plotpdftwo,plotlogden,plotpdftwom,0.1,0.9)
plotfittwom = sim_fitlnpdf(plotpdftwom,plotlogden,plotpdftwom,0.1,0.9)

xcap = 10**(0.0+0.05*6)
ycap = 10**(-6.0+0.9*4)
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
plt.plot(plotden,plotfitm,'k--',plotden,plotfittwom,'r--')
pa, pb, = plt.plot(plotden,plotpdfm,'k',plotden,plotpdftwom,'r')
plt.legend((pa, pb), (r"$\displaystyle t = 0.59~t_{\rm ff}$",r"$\displaystyle t = 1.06~t_{\rm ff}$"),prop={'size':8})
plt.xscale('log')
plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0,1.0e6,1.0e-4,1.0e-2])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6])

plt.xlabel(r"$\displaystyle n_H / {\rm cm^{-3}}$")
plt.ylabel(r"$\displaystyle P_M$")

pp = PdfPages(plotdir + 'f5a.pdf')
pp.savefig()
pp.close()


