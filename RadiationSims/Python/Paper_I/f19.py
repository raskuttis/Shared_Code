from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

mlist = ['5.0e3', '1.0e4', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '15.0', '20.0', '20.0', '25.0', '35.0', '15.0', '25.0', '35.0']
nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflistnf = np.core.defchararray.add(dflist, ['_N256_Tf4_NF/'] * nfs)
dflist = np.core.defchararray.add(dflist, ['_N256_Tf4/'] * nfs)

modelsigma = []
modeleff = []
modeleffof = []

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    datafoldernf = dflistnf[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    hstnfdata = read_hstfile(hostname,datadir + datafoldernf + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    
    eff = hst_eff(hstdata, mcloud)
    eof = max(hst_mof(hstnfdata, mcloud, tff))
    eof = 0.12
    
    modeleffof.append(eff / (1.0 - eof))
    modeleff.append(eff)
    modelsigma.append(sigmacloud)

meff,beff = np.polyfit(np.log10(modelsigma), modeleff, 1)
mof,bof = np.polyfit(np.log10(modelsigma), modeleffof, 1)

alist = ['_a0.1_RT/', '_a0.2_RT/', '_a0.4_RT/', '_a0.8_RT/', '_a1.5_RT/', '/', '_a3.0_RT/', '_a6.0_RT/', '_a10.0_RT/']
anflist = ['_a0.1_NF/', '_a0.2_NF/', '_a0.4_NF/', '_a0.8_NF/', '_a1.5_NF/', '_NF/', '_a3.0_NF/', '_a6.0_NF/', '_a10.0_NF/']
nfs = len(alist)
dflist = np.core.defchararray.add(['UV_M5.0e4_R15.0_N256_Tf4'] * nfs, alist)
dflistnf = np.core.defchararray.add(['UV_M5.0e4_R15.0_N256_Tf4'] * nfs, anflist)

modelalpha = []
amodeleff = []
amodeleffof = []

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    datafoldernf = dflistnf[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    hstnfdata = read_hstfile(hostname,datadir + datafoldernf + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    alphavir = out_alphavir(outlines)
    sigmacloud = out_sigma(outlines)
    
    eff = hst_eff(hstdata, mcloud)
    eof = max(hst_mof(hstnfdata, mcloud, tff))
    
    amodeleffof.append(eff / (1.0 - eof))
    amodeleff.append(eff)
    modelalpha.append(alphavir)

sigmafit = np.logspace(1.0, 3.0, num=100, base=10.0)
efffit = meff*np.log10(sigmafit)+beff
eoffit = mof*np.log10(sigmafit)+bof

meff,beff = np.polyfit(np.log10(modelalpha), amodeleff, 1)
mof,bof = np.polyfit(np.log10(modelalpha), amodeleffof, 1)
alphafit = np.logspace(-1.0, 1.0, num=100, base=10.0)
aefffit = meff*np.log10(alphafit)+beff
aeoffit = mof*np.log10(alphafit)+bof

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plt.plot(modelsigma,modeleff,'k.',modelsigma,modeleffof,'r.',sigmafit,efffit,'k',sigmafit,eoffit,'r')
plt.text(10**(1+0.05*2),0.9,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1])
plt.yticks([0,0.5,1])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle \varepsilon$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(1,2,2)
plt.plot(modelalpha,amodeleff,'k.',modelalpha,amodeleffof,'r.',alphafit,aefffit,'k',alphafit,aeoffit,'r')
plt.text(10**(-1+0.05*2),0.9,r"$\displaystyle(b)$")
plt.xscale('log')
plt.axis([1.0e-1,1.0e1,0,1])
plt.yticks([0,0.5,1], [' ', ' ', ' '])
plt.xticks([1.0e-1,1.0,1.0e1])

#plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")
plt.xlabel(r"$\displaystyle \alpha_{\rm vir,0}$")

pp = PdfPages(plotdir + 'f19.pdf')
pp.savefig()
pp.close()


