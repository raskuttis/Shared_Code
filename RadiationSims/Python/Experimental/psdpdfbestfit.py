from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_pdf import *
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

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs0.4/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs0.5/',
          'UV_M5.0e4_R15.0_N256_Tf4_Cs0.75/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs1.0/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs1.5/',
          'UV_M5.0e4_R15.0_N256_Tf4_Cs2.0/']
nds = len(dflist)

plotmachs = []
plotsigmas = []
plotmeans = []
plotsurfds = []
plotmacherrors = []

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
    csound = out_csound(outlines)
    machn = hst_mach(hstdata, csound)
    surfd = out_sigma(outlines)
    nlines = len(pdflines)
    nts = int(np.floor(nlines / 2))
    times = []
    sigmas = []
    means = []
    pdfmachs = []
    
    eff = hst_eff(hstdata, mcloud)
    tmini = np.argmax(mstar > 0.1 * eff)
    tmin = np.max([time[tmini],0.5])
    tmin = tmin/tmin*0.49
    tminpdfi = np.abs(pdftime - tmin).argmin()
    tmaxi = np.argmax(mstar > 0.9 * eff)
    tmax = np.max([time[tmaxi],1.5])
    tmax = tmax/tmax*0.51
    tmaxpdfi = np.abs(pdftime - tmax).argmin()
    
    for j in xrange(tminpdfi, tmaxpdfi):
        
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
        pdfmachs.append(machn[np.argmax(time > pdftime[j])])
        times.append(pdftime[j])
        sigmas.append(sigmam)
        means.append(meanm)
    
    sigmean = np.mean(sigmas)
    sigstd = np.std(sigmas)
    machmean = np.mean(pdfmachs)
    machsigma = np.std(pdfmachs)
    print machmean, machsigma, sigmean, sigstd

    plotsigmas.append(sigstd)
    plotmeans.append(sigmean)
    plotmachs.append(machmean)
    plotmacherrors.append(machsigma)
    plotsurfds.append(surfd)

modelmachn = np.logspace(0.0,2.0,num=100)
modelsigma = np.sqrt(np.log(1.0 + modelmachn**2/4.0))

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
plt.errorbar(plotmachs,plotmeans,xerr=[plotmacherrors,plotmacherrors],yerr=[plotsigmas,plotsigmas],fmt='.')
plt.plot(modelmachn,modelsigma,'k--')
#plt.text(0.05*3,0.9*3,r"$\displaystyle(a)$")
plt.axis([1.0,30.0,1,3])
plt.xticks([1.0,1.0e1])
plt.yticks([0,1,2,3])
plt.xscale('log')

plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")
plt.xlabel(r"$\displaystyle \mathcal{M}$")

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()


