from matplotlib.backends.backend_pdf import PdfPages
from hyp_math import *
from hyp_read import *
from hyp_out import *
from hyp_hst import *
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


modelsigma = np.logspace(1.0, 3.0, num = 1000.0)
loweff = []
mideff = []
higheff = []
fideff = []
x = 1.0
psi = 2000.0
nfs = len(modelsigma)
for i in xrange(0,nfs):
    sigma0 = modelsigma[i]
    sigma = 0.0
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    loweff.append(min_func(eps))
    sigma = 0.5
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    mideff.append(min_func(eps))
    sigma = 1.0
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    higheff.append(min_func(eps))
    sigma = 1.5
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    fideff.append(min_func(eps))

lowmsigma = np.logspace(1.0, 3.0, num = 15.0)
lowmeff = []
nfs = len(lowmsigma)
for i in xrange(0,nfs):
    sigma0 = lowmsigma[i]
    sigma = 0.0
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    lowmeff.append(min_func(eps))

psimodelsigma = np.logspace(1.0, 3.0, num = 1000.0)
psiloweff = []
psimideff = []
psihigheff = []
x = 1.0
psi = 2000.0
sigma = 1.5
nfs = len(modelsigma)
for i in xrange(0,nfs):
    sigma0 = modelsigma[i]
    psi = 400.0
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    psiloweff.append(min_func(eps))
    psi = 800.0
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    psimideff.append(min_func(eps))
    psi = 2000.0
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    psihigheff.append(min_func(eps))

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_P200/', 'UV_M5.0e4_R15.0_N256_Tf4_P800/',
          'UV_M5.0e4_R15.0_N256_Tf4/']
nds = len(dflist)
sigmas = []
effs = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    eff = hst_eff(hstdata, mcloud)
    sigmacl = out_sigma(outlines)
    effs.append(eff)
    sigmas.append(sigmacl)

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plow, pmid, phigh, pfid = plt.plot(modelsigma,loweff,'k',modelsigma,mideff,'r',modelsigma,higheff,'b',modelsigma,fideff,'g')
plt.plot(lowmsigma,lowmeff,'r.')
plt.text(10**(1+0.9*2),0.9,r"$\displaystyle(a)$")
plt.legend((plow, pmid, phigh, pfid), (r"$\displaystyle \sigma_{{\rm ln}\Sigma} = 0.0$",r"$\displaystyle 0.5$",r"$\displaystyle 1.0$",r"$\displaystyle 1.5$"),prop={'size':8},loc=2)
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1])
plt.yticks([0,0.5,1])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle \varepsilon$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(1,2,2)
plow, pmid, phigh, = plt.plot(psimodelsigma,psiloweff,'k',psimodelsigma,psimideff,'r',psimodelsigma,psihigheff,'b')
plt.plot(sigmas[0], effs[0], 'k.', sigmas[1], effs[1], 'r.', sigmas[2], effs[2], 'b.')
plt.text(10**(1+0.05*2),0.9,r"$\displaystyle(b)$")
plt.legend((plow, pmid, phigh), (r"$\displaystyle \Psi = 200.0$",r"$\displaystyle 800.0$",r"$\displaystyle 2000.0$"),prop={'size':8},loc=4)
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1])
plt.xticks([1.0e1,1.0e2,1.0e3])
plt.yticks([0,0.5,1], [' ', ' ', ' '])

#plt.ylabel(r"$\displaystyle \varepsilon$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'f28.pdf')
pp.savefig()
pp.close()


