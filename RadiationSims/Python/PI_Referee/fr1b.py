from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_star import *
from hyp_pdf import *
import numpy as np

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
pdffile = 'denpdfall.dat'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs0.5/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs1.0/',
          'UV_M5.0e4_R15.0_N256_Tf4_Cs2.0/']
nds = len(dflist)

plotdens = []
plotpdfs = []
plotpdfms = []
altplotpdfs = []
altplotpdfms = []

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

    pdftime = sim_pdftime(pdflines, tff)
    if (i == 0):
        eff = hst_eff(hstdata, mcloud)
        tmini = np.abs(mstar - 0.1 * eff).argmin()
        tmin = time[tmini]
        thmini = np.abs(mstar - 0.9 * eff).argmin()
        thmin = time[thmini]
    
    tminpdfi = np.abs(pdftime - tmin).argmin()
    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotpdfm = sim_pdfm(pdflines,tminpdfi)
    plotpdf = plotpdf / np.sum(plotpdf)

    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden

    plotpdfs.append(plotpdf)
    plotdens.append(plotden)
    plotpdfms.append(plotpdfm)

    tminpdfi = np.abs(pdftime - thmin).argmin()

    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotpdfm = sim_pdfm(pdflines,tminpdfi)
    plotpdf = plotpdf / np.sum(plotpdf)
    
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    
    altplotpdfs.append(plotpdf)
    altplotpdfms.append(plotpdfm)


xcap = 10**(0.05*6)
ycap = 10**(-8+0.9*6)
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
#plt.subplots_adjust(left=0.2)
for i in xrange(0,nds):
    plt.plot(plotdens[i],plotpdfs[i])
    plt.plot(plotdens[i],plotpdfms[i],'--')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0,1.0e6,1.0e-8,1.0e-2])
plt.yticks([1.0e-8,1.0e-6,1.0e-4,1.0e-2])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6])

plt.xlabel(r"$\displaystyle n_H / {\rm cm^{-3}}$")
plt.ylabel(r"$\displaystyle P$")

plt.subplot(1,2,2)
for i in xrange(0,nds):
    plt.plot(plotdens[i],altplotpdfs[i])
    plt.plot(plotdens[i],altplotpdfms[i],'--')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(b)$")
plt.axis([1.0,1.0e6,1.0e-8,1.0e-2])
plt.yticks([])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6])

plt.xlabel(r"$\displaystyle n_H / {\rm cm^{-3}}$")
#plt.ylabel(r"$\displaystyle P$")

pp = PdfPages(plotdir + 'fr1b.pdf')
pp.savefig()
pp.close()


