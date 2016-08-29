from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import fminbound

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

# Missing 1e5r15

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

x = 0.84
sigma = 1.42
sigmalow = 0.3
psi = 2000.0

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    datafoldernf = dflistnf[i]
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
    modelsigma.append(sigmacloud * (1.0 - 0.12) / x**2)
    print i, datafolder, eff, eof, eff / (1.0 - eof)


akt = 2.5
machkt = 30.0
effkt = 0.01
sigkt = ktvals(akt, machkt)
fitsigma = np.logspace(0.5, 3.5, num = 1000.0)
fiteff = []
fitlowsigeff = []
ktfiteff = []
ktmodelfiteff = []
nfs = len(fitsigma)
for i in xrange(0,nfs):
    sigma0 = fitsigma[i]
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    fiteff.append(min_func(eps))
    def min_func_lows(eps):
        return 1.0 - eps_of(eps, x, sigmalow, sigma0, psi)
    eps = fminbound(min_func_lows,1.0e-6,1.0)
    fitlowsigeff.append(min_func_lows(eps))
    ktfiteff.append(kteff(psi,effkt,sigma0,sigkt))
    ktmodelfiteff.append(kteff(psi,0.44,sigma0,sigma))

nconv = 50
ktmodelfiteff = np.convolve(ktmodelfiteff, np.ones((nconv,))/nconv, mode='same')
sfiteff = eff_simple(1.0, 2000.0, fitsigma)

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.plot(modelsigma,modeleffof,'r.',fitsigma,fiteff,'r',fitsigma,sfiteff,'k',fitsigma,fitlowsigeff,'--r',)
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1])
plt.yticks([0,0.2,0.4,0.6,0.8,1])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle \varepsilon$")
plt.xlabel(r"$\displaystyle \langle \Sigma_{\rm cloud} \rangle / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'fpresd.pdf')
pp.savefig()
pp.close()




