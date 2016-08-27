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
pdffile = 'sdpdfxy.dat'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4/']
pdflist = ['sdpdfxy.dat', 'sdpdfmeancirc.dat']
nds = len(dflist)

plottimes = []
plotsigmas = []
plotmeans = []
plotxs = []

for i in xrange(0,nds):
    
    datafolder = dflist[i]
    pdffile = pdflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    msol = out_msol(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    pdftime = sim_pdftime(pdflines, tff)
    nlines = len(pdflines)
    nts = int(np.floor(nlines / 2))
    nts = 5000
    times = []
    sigmas = []
    means = []
    
    for j in xrange(300, nts-1):
        
        pdf = sim_pdf(pdflines,j)
        nvals = np.sum(pdf)
        pdflogx = sim_pdfx(pdflines)
        pdfm = sim_pdfm(pdflines,j)
        #amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.1,0.9)
        ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.1,0.9,nvals)
        #sigma = sigma / np.log10(np.exp(1))
        #mu = mu / np.log10(np.exp(1))
        sigmam = sigmam / np.log10(np.exp(1))
        mum = mum / np.log10(np.exp(1))
        #mean = np.exp(mu + 0.5 * sigma**2)
        meanm = np.exp(mum - 0.5 * sigmam**2)
        
        print j, pdftime[j], mum, sigmam, meanm
        
        if sigmam > 0.0:
            times.append(pdftime[j])
            sigmas.append(sigmam)
            means.append(meanm / msol)

    plottimes.append(times)
    plotsigmas.append(sigmas)
    plotmeans.append(means)


plt.figure(figsize = [xsize,3*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(3,1,1)
plt.subplots_adjust(left=0.2)
pfid, pgrav, = plt.plot(plottimes[0],plotsigmas[0],'k',plottimes[1],plotsigmas[1],'r')
plt.legend((pfid, pgrav), (r"Fiducial",r"Circumcluster"))
plt.text(0.05*3,0.9*3,r"$\displaystyle(a)$")
plt.axis([0,3,0,3])
plt.yticks([0,1,2,3])
#plt.xticks([0,1,2,3])
plt.xticks([])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")

plt.subplot(2,1,2)
#plt.subplots_adjust(left=0.2)
plt.plot(plottimes[0],plotmeans[0],'k',plottimes[1],plotmeans[1],'r')
plt.text(0.05*3,0.9*10,r"$\displaystyle(b)$")
plt.axis([0,3,0.0,60.0])
plt.yticks([0.0,10.0,20.0,30.0,40.0,50.0,60.0])
plt.xticks([0,1,2,3])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \langle \Sigma \rangle_{\rm cloud}$")

pp = PdfPages(plotdir + 'fr3circ.pdf')
pp.savefig()
pp.close()


