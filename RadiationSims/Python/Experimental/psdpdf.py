from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_pdf import *
import numpy as np

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
pdffile = 'sdpdfxy.dat'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N128_Tf4_B0.1_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_Ds0.5/', 'UV_M5.0e4_R15.0_N256_Tf4_Ds0.8/', 'UV_M5.0e4_R15.0_N256_Tf4_Ds1.0/']
nds = len(dflist)

plotdens = []
plotpdfs = []
plotpdfms = []
nconv = 25

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
    msol = out_msol(outlines)

    eff = hst_eff(hstdata, mcloud)
    pdftime = sim_pdftime(pdflines, tff)

    tmini = np.abs(time - [1.2]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotpdfm = sim_pdfm(pdflines,tminpdfi)
    plotpdf = plotpdf / np.sum(plotpdf)
    plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
    plotpdfm = np.convolve(plotpdfm, np.ones((nconv,))/nconv, mode='same')
    #print plotpdf

    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    #print plotden

    plotpdfs.append(plotpdf)
    plotdens.append(plotden)
    plotpdfms.append(plotpdfm)


plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
colors = ['k','r','b','g','m','c']
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
for i in xrange(0,nds):
    plt.plot(plotdens[i],plotpdfms[i],colors[i])
#plt.plot(plotdens[i],plotpdfms[i],'-.')
plt.xscale('log')
plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])

plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_M$")

pp = PdfPages(plotdir + 'psdpdf.pdf')
pp.savefig()
pp.close()


