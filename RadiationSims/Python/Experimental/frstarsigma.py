from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_fluxes import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_pdf import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'sdpdfxy.dat'
starfile = 'star'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

mlist = ['5.0e3', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '15.0', '20.0', '15.0', '20.0', '25.0', '35.0', '15.0', '25.0', '35.0']
nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflist = np.core.defchararray.add(dflist, ['_N256_Tf4/'] * nfs)

modelsigma = []
modeleff = []
modeleffof = []
tplot = [0.1,0.5,0.9]
ntplot = len(tplot)
nconvpdf = 50

allfabs = np.zeros((3,nfs))
allsigma = np.zeros(nfs)
allfedd = np.zeros((3,nfs))

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    psi = out_Psi(outlines)
    clight = out_clight(outlines)
    kappa = out_kappa(outlines)
    rcloud = out_rcloud(outlines)
    sigmaadj = sigmacloud * 0.88
    
    time = hst_time(hstdata, tff)
    eff = hst_mstar(hstdata, mcloud)
    eps = hst_eff(hstdata, mcloud)
    print i, datafolder, max(time), sigmacloud
    
    stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)
    dcom = star_sigma(stardata) / rcloud
    
    allsigma[i] = sigmacloud
    
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    pdftime = sim_pdftime(pdflines, tff)
    nts = len(pdftime)
    ntmax = np.min([7500,nts-100])

    plottime = []
    plotx = []
    startflag = False
    for j in xrange(100,ntmax):
    
        epst = eff[j]
        plotpdf = sim_pdf(pdflines,j)
        plotpdf = np.convolve(plotpdf, np.ones((nconvpdf,))/nconvpdf, mode='same')
        plotpdfm = sim_pdfm(pdflines,j)
        plotpdfm = np.convolve(plotpdfm, np.ones((nconvpdf,))/nconvpdf, mode='same')
        plotpdf = plotpdf / np.sum(plotpdf)
        nvals = np.sum(plotpdf)
    
        plotlogden = sim_pdfx(pdflines)
        plotden = 10**plotlogden / msol
        if startflag:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdfm,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
            ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(plotpdfm,plotlogden,plotpdfm,0.15,0.95,nvals,pin=[ampmold,mumold,sigmamold])
        else:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdfm,0.15,0.95,nvals)
            ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(plotpdfm,plotlogden,plotpdfm,0.15,0.95,nvals)
    
        if (sigma > 0.0):
        
            [ampold, sigmaold, muold] = [amp, sigma, mu]
            [ampmold, sigmamold, mumold] = [ampm, sigmam, mum]
            sigma = sigma / np.log10(np.exp(1))
            mu = mu / np.log10(np.exp(1))
            sigmam = sigmam / np.log10(np.exp(1))
            mum = mum / np.log10(np.exp(1))
            mean = np.exp(mu + 0.5 * sigma**2)
            meanm = np.exp(mum - 0.5 * sigmam**2)
            x = np.sqrt(sigmaadj * msol * (1.0 - epst) / (mean * 1.0))
            xm = np.sqrt(sigmaadj * msol * (1.0 - epst) / (meanm * 1.0))
            startflag = True
        
            plottime.append(pdftime[j])
            plotx.append(x)

    for j in xrange(0,ntplot):
        
        tmini = np.abs(eff - tplot[j] * eps).argmin()
        thist = time[tmini]
        tminfi = np.abs(plottime - thist).argmin()
        dscom = star_masssigma(stardata,time,thist) / rcloud
        allfabs[j,i] = dscom
        allfedd[j,i] = plotx[tminfi]
    
        print j, dcom[tmini], allfedd[j,i], dscom

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plow, pmid, phigh = plt.plot(allsigma,allfabs[0,:],'k.',allsigma,allfabs[1,:],'r.',allsigma,allfabs[2,:],'b.')
plt.legend((plow,pmid,phigh), (r"$\displaystyle 10~\%$",r"$\displaystyle 50~\%$",r"$\displaystyle 90~\%$"),prop={'size':8})
plt.text(10**(1+0.05*2),0.9*2,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,2])
plt.yticks([0,0.5,1,1.5,2.0])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

axfedd = plt.subplot(1,2,2)
axfedd.yaxis.set_label_position("right")
plow, pmid, phigh = plt.plot(allsigma,allfedd[0,:],'k.',allsigma,allfedd[1,:],'r.',allsigma,allfedd[2,:],'b.')
plt.text(10**(1+0.05*2),0.9*2.0,r"$\displaystyle(b)$")
plt.xscale('log')
#plt.yscale('log')
plt.axis([1.0e1,1.0e3,0.0,2.0])
plt.yticks([0,0.5,1,1.5,2.0],[' ',' ',' ',' '])
plt.xticks([1.0e1,1.0e2,1.0e3])

#plt.ylabel(r"$\displaystyle f_{\rm edd,cum}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'fstarsigma.pdf')
pp.savefig()
pp.close()

