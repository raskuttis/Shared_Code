from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_math import *
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

mlist = ['5.0e3', '1.0e4', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '2.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '15.0', '20.0', '20.0', '25.0', '15.0', '25.0', '35.0']
nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflist = np.core.defchararray.add(dflist, ['_N256_Tf4/'] * nfs)

modelsigma = []
tstars = []
tbreaks = []
mbreaks = []
bplbetas = []
plbetas = []
bplchis = []
bplchilims = []

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    tstar = time[np.argmax(mstar > 0.0)]
    eff = hst_eff(hstdata, mcloud)
    if sigmacloud >= 52:
        maxm = min([0.3,0.9*eff])
    else:
        maxm = min([0.2,0.9*eff])
    
    tcut, afpl, aspl, anfpl, anspl = mstar_brokenplfit(time, mstar, 0.01, maxm)
    pfit, bfit = mstar_plfit(time, mstar, 0.01, maxm)
    mcut = mstar[np.argmax(time > tcut)]
    
    mstarfit = broken_pl(time, tstar, max(time)-0.01, pfit, 1.0, 10**bfit, 1.0)
    plmstarfit = broken_pl(time, tstar, tcut, afpl, aspl, anfpl, anspl)

    nsq = len(time[np.argmax(mstar > 0.01):np.argmax(mstar > maxm)])
    nsqlim = len(time[np.argmax(mstar > mcut):np.argmax(mstar > maxm)])
    chisq = lsqerror(mstar, plmstarfit, 0.01, maxm) / (lsqerror(mstar, mstarfit, 0.01, maxm))
    chisqlim = (lsqerror(mstar, plmstarfit, mcut, maxm) * nsq) / (lsqerror(mstar, mstarfit, mcut, maxm) * nsqlim)

    print i, datafolder, tcut, aspl

    modelsigma.append(sigmacloud)
    tstars.append(tstar)
    tbreaks.append(tcut)
    plbetas.append(pfit)
    bplbetas.append(aspl)
    bplchis.append(chisq)
    bplchilims.append(chisqlim)
    mbreaks.append(mcut)

plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.plot(modelsigma,tstars,'k.',modelsigma,tbreaks,'r.')
plt.text(10**(1+0.05*2),0.9*1.5,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,1.5])
plt.yticks([0,0.5,1,1.5])
#plt.xticks([1.0e1,1.0e2,1.0e3])
plt.xticks([])

plt.ylabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

axbeta = plt.subplot(2,2,2)
plt.plot(modelsigma,plbetas,'k.',modelsigma,bplbetas,'r.')
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axbeta.yaxis.set_label_position("right")
plt.text(10**(1+0.05*2),0.9*1.5+0.5,r"$\displaystyle(b)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0.5,2.0])
plt.yticks([0.5,1,1.5,2.0])
#plt.xticks([1.0e1,1.0e2,1.0e3])
plt.xticks([])

plt.ylabel(r"$\displaystyle \beta$")
#plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(2,2,3)
plt.plot(modelsigma,mbreaks,'k.')
plt.text(10**(1+0.05*2),0.9*0.3,r"$\displaystyle(c)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,0.3])
plt.yticks([0,0.1,0.2,0.3])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

axchi = plt.subplot(2,2,4)
plt.plot(modelsigma,bplchis,'k.',modelsigma,bplchilims,'r.')
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axchi.yaxis.set_label_position("right")
plt.text(10**(1+0.05*2),10**(-2+0.9*2),r"$\displaystyle(d)$")
plt.xscale('log')
plt.yscale('log')
plt.axis([1.0e1,1.0e3,1.0e-2,1.0])
plt.yticks([1.0e-2,1.0e-1,1.0])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle \chi^2_{\rm broken} / \chi^2_{\rm single}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'f24.pdf')
pp.savefig()
pp.close()


