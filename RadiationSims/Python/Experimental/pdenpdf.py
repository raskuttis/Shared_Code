from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_pdf import *
import numpy as np

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
pdffile = 'denpdfall.dat'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B0.2_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.02_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.002_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.0002_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.00002_SDs/']
nds = len(dflist)

plotdens = []
plotpdfs = []
plotpdfms = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)

    eff = hst_eff(hstdata, mcloud)
    pdftime = sim_pdftime(pdflines, tff)

    tmini = np.abs(time - [0.6]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotpdfm = sim_pdfm(pdflines,tminpdfi)
    plotpdf = plotpdf / np.sum(plotpdf)

    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden

    plotpdfs.append(plotpdf)
    plotdens.append(plotden)
    plotpdfms.append(plotpdfm)


plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

colors = ['k','r','b','g','c']
plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
for i in xrange(0,nds):
    plt.plot(plotdens[i],plotpdfms[i],colors[i])
#plt.plot(plotdens[i],plotpdfms[i],colors[i],'-.')
plt.xscale('log')
plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0,1.0e6,1.0e-8,1.0e-2])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6])

plt.xlabel(r"$\displaystyle n_H / {\rm cm^{-3}}$")
plt.ylabel(r"$\displaystyle P$")

pp = PdfPages(plotdir + 'pdenpdf.pdf')
pp.savefig()
pp.close()


