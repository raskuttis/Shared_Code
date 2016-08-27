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

dflist = ['UV_M5.0e4_R15.0_N128_Tf4_Ds0.1/', 'UV_M5.0e4_R15.0_N128_Tf4_Ds0.3/', 'UV_M5.0e4_R15.0_N128_Tf4_Ds0.5/',
          'UV_M5.0e4_R15.0_N128_Tf4_Ds0.8/', 'UV_M5.0e4_R15.0_N128_Tf4_Ds1.0/']
nds = len(dflist)

plottimes = []
plotsigmas = []
plotmeans = []

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
    nlines = len(pdflines)
    nts = int(np.floor(nlines / 2))
    times = []
    sigmas = []
    means = []
    
    for j in xrange(300, 6000):
        
        pdf = sim_pdf(pdflines,j)
        pdflogx = sim_pdfx(pdflines)
        pdfm = sim_pdfm(pdflines,j)
        #amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.1,0.9)
        ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.1,0.9)
        #sigma = sigma / np.log10(np.exp(1))
        #mu = mu / np.log10(np.exp(1))
        sigmam = sigmam / np.log10(np.exp(1))
        mum = mum / np.log10(np.exp(1))
        #mean = np.exp(mu + 0.5 * sigma**2)
        meanm = np.exp(mum - 0.5 * sigmam**2)
        times.append(pdftime[j])
        sigmas.append(sigmam)
        means.append(meanm)

    plottimes.append(times)
    plotsigmas.append(sigmas)
    plotmeans.append(means)


plt.figure(figsize = [xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

colors = ['k','r','b','g','c']
plt.subplot(2,1,1)
plt.subplots_adjust(left=0.2)
for i in xrange(0,nds):
    plt.plot(plottimes[i],plotsigmas[i],colors[i])
plt.text(0.05*3,0.9*3,r"$\displaystyle(a)$")
plt.axis([0,3,0,3])
plt.yticks([0,1,2,3])
#plt.xticks([0,1,2,3])
plt.xticks([])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\rho}$")

plt.subplot(2,1,2)
#plt.subplots_adjust(left=0.2)
for i in xrange(0,nds):
    plt.plot(plottimes[i],plotmeans[i],colors[i])
plt.text(0.05*3,10**(0.9*3),r"$\displaystyle(b)$")
plt.axis([0,3,1.0,1.0e3])
plt.yscale('log')
plt.yticks([1.0,1.0e1,1.0e2,1.0e3])
plt.xticks([0,1,2,3])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \langle \rho \rangle$")

pp = PdfPages(plotdir + 'fr3c.pdf')
pp.savefig()
pp.close()


