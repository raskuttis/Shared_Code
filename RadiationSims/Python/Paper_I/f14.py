from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy import integrate
from scipy import special

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

# Missing 5e4r15 and 2e5r25

mlist = ['5.0e3', '1.0e4', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '20.0', '15.0', '20.0', '25.0', '35.0', '15.0', '35.0']
nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflistnf = np.core.defchararray.add(dflist, ['_N256_Tf4_NF/'] * nfs)

modelsigma = []
modelof = []

for i in xrange(0,nfs):
    
    datafoldernf = dflistnf[i]
    outlines = read_outfile(hostname,datadir + datafoldernf + outfile)
    hstdata = read_hstfile(hostname,datadir + datafoldernf + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    
    time = hst_time(hstdata, tff)
    mofall = hst_mofsimple(hstdata, mcloud)
    mofold = max(mofall)
    mofall = hst_mof(hstdata, mcloud, tff)
    mof = max(mofall)
    print mof, mofold, tff, mcloud
    
    modelof.append(mof)
    modelsigma.append(sigmacloud)
    print i, datafoldernf, mof, max(time)

anflist = ['_a0.1_NF/', '_a0.2_NF/', '_a0.4_NF/', '_a0.8_NF/', '_a1.5_NF/', '_NF/', '_a3.0_NF/', '_a6.0_NF/', '_a10.0_NF/']
nfs = len(anflist)
dflistnf = np.core.defchararray.add(['UV_M5.0e4_R15.0_N256_Tf4'] * nfs, anflist)

modelalpha = []
modelaof = []

for i in xrange(0,nfs):
    
    datafoldernf = dflistnf[i]
    outlines = read_outfile(hostname,datadir + datafoldernf + outfile)
    hstdata = read_hstfile(hostname,datadir + datafoldernf + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    alphavir = out_alphavir(outlines)
    
    time = hst_time(hstdata, tff)
    mofall = hst_mofsimple(hstdata, mcloud)
    mofold = max(mofall)
    mofall = hst_mof(hstdata, mcloud, tff)
    mof = max(mofall)
    print mof, mofold
    
    modelaof.append(mof)
    modelalpha.append(alphavir)
    print i, datafoldernf, mof

alphafit = np.logspace(-1.0, 1.0, num=100, base=10.0)
def integrand(x, alpha):
    return 3.0 * x**2 * (1.0 - special.erf(np.sqrt(5.0*(3-x**2)/(6*alpha))))
efffit = []
for i in xrange(0,100):
    alpha = alphafit[i]
    eff = integrate.quad(integrand, 0, 1, args = alpha)
    efffit.append(eff[0])
def integrand(x, alpha):
    return 3.0 * x**2 * (1.0 - special.erf(np.sqrt(5.0*(2-x**2)/(6*alpha))))
twoefffit = []
for i in xrange(0,100):
    alpha = alphafit[i]
    eff = integrate.quad(integrand, 0, 1, args = alpha)
    twoefffit.append(eff[0])


plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plt.plot(modelsigma,modelof,'r.')
plt.text(10**(1+0.05*2),10**(-3+0.9*3),r"$\displaystyle(a)$")
plt.xscale('log')
plt.yscale('log')
plt.axis([1.0e1,1.0e3,1.0e-3,1])
plt.yticks([1.0e-3,1.0e-2,1.0e-1,1.0])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle \varepsilon_{of}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(1,2,2)
plt.plot(modelalpha,modelaof,'r.',alphafit,efffit,'k',alphafit,twoefffit,'r')
plt.text(10**(-1+0.05*2),10**(-3+0.9*3),r"$\displaystyle(b)$")
plt.xscale('log')
plt.yscale('log')
plt.axis([1.0e-1,1.0e1,1.0e-3,1.0])
plt.yticks([])
plt.xticks([1.0e-1,1.0,1.0e1])

#plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")
plt.xlabel(r"$\displaystyle \alpha_{\rm vir,0}$")

pp = PdfPages(plotdir + 'f17.pdf')
pp.savefig()
pp.close()


