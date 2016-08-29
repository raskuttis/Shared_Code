from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_pdf import *
import numpy as np

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
pdffile = 'sdmasspdfmeancirc.dat'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

pdflist = ['sdpdfxy.dat', 'sdpdfyz.dat', 'sdmasspdfallcirc.dat', 'sdmasspdfmeancirc.dat']
pdftype = [1,1,0,0]
nps = len(pdflist)
dflist = ['UV_M5.0e4_R15.0_N128_Tf4_NF_B0.05_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_NF_B0.5_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_NF_B5.0_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_NF_B50.0_SDs/']
betas = [0.05, 0.5, 5.0, 50.0]
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
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    msol = out_msol(outlines)
    gravc = out_G(outlines)
    
    #bmag = np.sqrt(2.0 * 0.2**2 / betas[i])
    #phimag = np.pi * rcloud**2 * bmag
    #mmag = phimag / (2.0 * np.pi * np.sqrt(gravc))
    #mtflux = mcloud / mmag
    #print mtflux

    eff = hst_eff(hstdata, mcloud)
    
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    pdftime = sim_pdftime(pdflines, tff)

    tmini = np.abs(time - [1.0]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotpdfm = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / np.sum(plotpdf)
    plotpdfm = plotpdfm / np.sum(plotpdfm)
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
pa, pb, pc, pd = plt.plot(plotdens[0],plotpdfms[0],colors[0],plotdens[1],plotpdfms[1],colors[1],plotdens[2],plotpdfms[2],colors[2],plotdens[3],plotpdfms[3],colors[3])
plt.legend((pa, pb, pc, pd), (r"$\displaystyle M / \Phi = 55$",r"$17$",r"$5.5$",r"$2$"),prop={'size':8},loc=2)
#plt.plot(plotdens[i],plotpdfms[i],'-.')
plt.xscale('log')
plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])

plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_M$")

pp = PdfPages(plotdir + 'psdpdfmag.pdf')
pp.savefig()
pp.close()


