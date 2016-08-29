from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_math import *
from ..Hyperion.hyp_pdf import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'sdpdfxy.dat'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

mlist = ['5.0e3', '1.0e4', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '15.0', '20.0', '20.0', '25.0', '35.0', '15.0', '25.0']
nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflistnf = np.core.defchararray.add(dflist, ['_N256_Tf4_NF/'] * nfs)
dflist = np.core.defchararray.add(dflist, ['_N256_Tf4/'] * nfs)

modelsigma = []
modeleff = []
modeleffnf = []
modeltcuts = []
modelslins = []

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    datafoldernf = dflistnf[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    hstnfdata = read_hstfile(hostname,datadir + datafoldernf + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    rhocloud = out_rhocloud(outlines)
    sigmacloud = out_sigma(outlines)
    
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    eps = hst_eff(hstdata, mcloud)
    if sigmacloud >= 52:
        maxm = min([0.3,0.9*eps])
    else:
        maxm = min([0.2,0.9*eps])

    tcutlin, anflin, anslin = mstar_brokenlinfit(time, mstar, 0.01, maxm)
    if (tcutlin > 2.0):
        tcutlin = time[np.argmax(mstar > 0.1)]
    modeltcuts.append(tcutlin)
    modelslins.append(anslin)

    pdftime = sim_pdftime(pdflines, tff)
    tminpdfi = np.abs(pdftime - tcutlin).argmin()
    pdf = sim_pdf(pdflines,tminpdfi)
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,tminpdfi)
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.1,0.9,nvals)
    sigma = sigma / np.log10(np.exp(1))
    mu = mu / np.log10(np.exp(1))
    mean = np.exp(mu - 0.5 * sigma**2)
    mean = rhocloud
    eff = anslin * np.sqrt(rhocloud / mean)
    print i, datafolder, tcutlin, anslin, rhocloud, mean

    time = hst_time(hstnfdata, tff)
    mstar = hst_mstar(hstnfdata, mcloud)
    eps = hst_eff(hstnfdata, mcloud)
    maxm = min([0.3,0.9*eps])
    tcutlin, anflin, anslin = mstar_brokenlinfit(time, mstar, 0.01, maxm)

    mean = rhocloud
    effnf = anslin * np.sqrt(rhocloud / mean)

    modeleffnf.append(effnf)
    modeleff.append(eff)
    modelsigma.append(sigmacloud)

alist = ['_a0.1_RT/', '_a0.2_RT/', '_a0.4_RT/', '_a0.8_RT/', '_a1.5_RT/', '_a3.0_RT/', '_a6.0_RT/']
anflist = ['_a0.1_NF/', '_a0.2_NF/', '_a0.4_NF/', '_a0.8_NF/', '_a1.5_NF/', '_a3.0_NF/', '_a6.0_NF/']
nfs = len(alist)
dflist = np.core.defchararray.add(['UV_M5.0e4_R15.0_N256_Tf4'] * nfs, alist)
dflistnf = np.core.defchararray.add(['UV_M5.0e4_R15.0_N256_Tf4'] * nfs, anflist)

modelalpha = []
amodeleff = []
amodeleffnf = []
nconv = 25

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    datafoldernf = dflistnf[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    hstnfdata = read_hstfile(hostname,datadir + datafoldernf + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    alphavir = out_alphavir(outlines)
    sigmacloud = out_sigma(outlines)
    
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    eps = hst_eff(hstdata, mcloud)
    maxm = min([0.3,0.9*eps])
    tcutlin, anflin, anslin = mstar_brokenlinfit(time, mstar, 0.01, maxm)
    if (tcutlin > 2.0):
        tcutlin = time[np.argmax(mstar > 0.1)]
    #print tcutlin, anslin
    tmini = np.abs(time - tcutlin).argmin()
    epst = mstar[tmini]

    pdftime = sim_pdftime(pdflines, tff)
    tminpdfi = np.abs(pdftime - tcutlin).argmin()
    pdf = sim_pdf(pdflines,tminpdfi)
    pdf = np.convolve(pdf, np.ones((nconv,))/nconv, mode='same')
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,tminpdfi)
    pdfm = np.convolve(pdfm, np.ones((nconv,))/nconv, mode='same')
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals)
    sigma = sigma / np.log10(np.exp(1))
    mu = mu / np.log10(np.exp(1))
    mean = np.exp(mu - 0.5 * sigma**2)
    x = np.sqrt(sigmacloud * msol * (1.0 - epst) / mean) / 0.84
    meancloud = rhocloud / x**3
    eff = anslin * np.sqrt(rhocloud / meancloud)
    #print mean / msol, sigmacloud, epst, x
    
    time = hst_time(hstnfdata, tff)
    mstar = hst_mstar(hstnfdata, mcloud)
    eps = hst_eff(hstnfdata, mcloud)
    maxm = min([0.3,0.9*eps])
    tcutlin, anflin, anslin = mstar_brokenlinfit(time, mstar, 0.01, maxm)
    effnf = anslin * np.sqrt(rhocloud / meancloud)
    #print tcutlin, anslin
    print i, datafolder, eff, effnf
    
    amodeleffnf.append(effnf)
    amodeleff.append(eff)
    modelalpha.append(alphavir)

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plt.plot(modelsigma,modeleff,'k.',modelsigma,modeleffnf,'r.')
plt.text(10**(1+0.05*2),0.9,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1])
plt.yticks([0,0.25,0.5,0.75,1])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle \varepsilon_{\rm ff}(\langle \rho \rangle)$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(1,2,2)
plt.plot(modelalpha,amodeleff,'k.',modelalpha,amodeleffnf,'r.')
plt.text(10**(-1+0.05*2),0.9,r"$\displaystyle(b)$")
plt.xscale('log')
plt.axis([1.0e-1,1.0e1,0,1])
plt.yticks([0,0.25,0.5,0.75,1],[' ',' ',' ',' ',' '])
plt.xticks([1.0e-1,1.0,1.0e1])

plt.xlabel(r"$\displaystyle \alpha_{\rm vir,0}$")

pp = PdfPages(plotdir + 'f34.pdf')
pp.savefig()
pp.close()


