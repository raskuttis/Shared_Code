from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_star import *
from hyp_pdf import *
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

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B0.02_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.002_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.0002_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.00002_SDs/']
nds = len(dflist)

plottimes = []
plotsigmas = []
plotmeans = []
plotmachs = []
nconv = 25
nconvall = 3

for i in xrange(0,nds):
    
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    msol = out_msol(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    csound = out_csound(outlines)
    machn = hst_mach(hstdata, csound)
    pdftime = sim_pdftime(pdflines, tff)
    nlines = len(pdflines)
    nts = int(np.floor(nlines / 2))
    times = []
    sigmas = []
    means = []
    pdfmachs = []
    jmin = 30
    jmax = np.min([1000,np.max(nts)-60])
    print jmin,jmax
    
    plottime = []
    plotsigma = []
    plotmean = []
    
    startflag = False
    
    for j in xrange(jmin,jmax):
        
        plotpdf = sim_pdfm(pdflines,j)
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        plotpdf = plotpdf / np.sum(plotpdf)
        nvals = np.sum(plotpdf)

        plotlogden = sim_pdfx(pdflines)
        plotden = 10**plotlogden / msol
        if startflag:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
        else:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals)
        
        if (sigma > 0.0):
            
            [ampold, sigmaold, muold] = [amp, sigma, mu]
            sigma = sigma / np.log10(np.exp(1))
            mu = mu / np.log10(np.exp(1))
            mean = np.exp(mu - 0.5 * sigma**2)
            print j, pdftime[j], sigma
            
            plottime.append(pdftime[j])
            plotsigma.append(sigma)
            plotmean.append(mean / msol)
            
            startflag = True

    plottimes.append(plottime)
    plotsigmas.append(np.convolve(plotsigma, np.ones((nconvall,))/nconvall, mode='same'))
    plotmeans.append(np.convolve(plotmean, np.ones((nconvall,))/nconvall, mode='same'))

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

colors = ['k','r','b','g','c','m','y']
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
for i in xrange(0,nds):
    plt.plot(plottimes[i],plotsigmas[i],colors[i])
plt.text(0.05*3,0.9*3,r"$\displaystyle(a)$")
plt.axis([0,2,0,3])
plt.yticks([0,1,2,3])
plt.xticks([0,1,2])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")

pp = PdfPages(plotdir + 'psdpdftimemag.pdf')
pp.savefig()
pp.close()


