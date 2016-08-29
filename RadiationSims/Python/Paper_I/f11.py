from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_math import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

mlist = ['5.0e3', '5.0e3', '5.0e3', '5.0e3', '1.0e4', '1.0e4', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '8.0', '10.0', '15.0', '5.0', '8.0', '10.0', '15.0', '5.0', '8.0', '10.0', '15.0', '20.0', '25.0', '8.0', '10.0', '15.0', '20.0', '15.0', '20.0', '25.0', '35.0', '15.0', '20.0', '35.0']
nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflist = np.core.defchararray.add(dflist, ['_N128_Tf4_a0.1/'] * nfs)

modelsigma = []
modeleff = []

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)

    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    mgas = hst_num(hstdata, 2)
    eff = hst_eff(hstdata, mgas[0])
    
    modelsigma.append(sigmacloud)
    modeleff.append(eff)
    print sigmacloud, eff

psi = 2000.0
xfit = 1.0
def min_func(sigmadata, xf):
    eps = eff_simple(xf, psi, sigmadata)
    return eps

xfit, pcov = curve_fit(min_func, modelsigma, modeleff)
sigmafit = np.logspace(1.0, 3.0, num=100, base=10.0)
efffit = eff_simple(xfit, psi, sigmafit)
print xfit

xcap = 0.05*3
ycap = 0.9*3
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.plot(modelsigma,modeleff,'r.',sigmafit,efffit,'k')
#plt.text(10**(1+0.05*2),ycap,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1])
plt.yticks([0,0.5,1])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle \varepsilon$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'f11.pdf')
pp.savefig()
pp.close()


