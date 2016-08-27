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

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4_NF/']
nds = len(dflist)

plotdens = []
plotpdfs = []
plotpdfms = []
altplotpdfs = []
altplotpdfms = []
plotfitms = []
altplotfitms = []

for i in xrange(0,nds):
    nconv = 25
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)

    if (i == 0):
        eff = hst_eff(hstdata, mcloud)
    pdftime = sim_pdftime(pdflines, tff)

    tmini = np.abs(mstar - 0.5 * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    print tmin

    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotpdfm = sim_pdfm(pdflines,tminpdfi)
    plotpdf = plotpdf / np.sum(plotpdf)

    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    plotfitm = sim_fitlnpdf(plotpdfm,plotlogden,plotpdfm,0.1,0.9)

    plotpdfs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same'))
    plotdens.append(plotden)
    plotpdfms.append(np.convolve(plotpdfm, np.ones((nconv,))/nconv, mode='same'))
    plotfitms.append(plotfitm)

    tmini = np.abs(mstar - 0.9 * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    print tmin
    
    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotpdfm = sim_pdfm(pdflines,tminpdfi)
    plotpdf = plotpdf / np.sum(plotpdf)
    
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    plotfitm = sim_fitlnpdf(plotpdfm,plotlogden,plotpdfm,0.1,0.9)
    
    altplotpdfs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same'))
    altplotpdfms.append(np.convolve(plotpdfm, np.ones((nconv,))/nconv, mode='same'))
    altplotfitms.append(plotfitm)

xcap = 10**(0.05*6)
ycap = 10**(-4+0.9*2)
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
#plt.subplots_adjust(left=0.2)
#plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[1],plotpdfs[1],'r')
plt.plot(plotdens[0],plotpdfms[0],'k',plotdens[1],plotpdfms[1],'r')
plt.plot(plotdens[0],plotfitms[0],'k--',plotdens[1],plotfitms[1],'r--')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0,1.0e6,1.0e-4,1.0e-2])
plt.yticks([1.0e-4,1.0e-3,1.0e-2])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6])

plt.xlabel(r"$\displaystyle n_H / {\rm cm^{-3}}$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(1,2,2)
#plt.plot(plotdens[0],altplotpdfs[0],'k',plotdens[1],altplotpdfs[1],'r')
plt.plot(plotdens[0],altplotpdfms[0],'k',plotdens[1],altplotpdfms[1],'r')
plt.plot(plotdens[0],altplotfitms[0],'k--',plotdens[1],altplotfitms[1],'r--')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(b)$")
plt.axis([1.0,1.0e6,1.0e-4,1.0e-2])
plt.yticks([])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6])

plt.xlabel(r"$\displaystyle n_H / {\rm cm^{-3}}$")
#plt.ylabel(r"$\displaystyle P$")

pp = PdfPages(plotdir + 'fr2a.pdf')
pp.savefig()
pp.close()


