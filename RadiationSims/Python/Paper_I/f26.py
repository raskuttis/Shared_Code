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

mlist = ['5.0e3', '1.0e4', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '15.0', '20.0', '20.0', '25.0', '35.0', '15.0', '25.0', '35.0']
nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflist = np.core.defchararray.add(dflist, ['_N256_Tf4/'] * nfs)

modelsigma = np.zeros(nfs)
tstars = np.zeros(nfs)
tsmax = np.zeros(nfs)
tsplus = np.zeros(nfs)
tsminus = np.zeros(nfs)
tamax = np.zeros(nfs)
taplus = np.zeros(nfs)
taminus = np.zeros(nfs)
tconvMyr = np.zeros(nfs)
effs = np.zeros(nfs)

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    
    tMyr = out_tMyr(outlines)
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    tstar = time[np.argmax(mstar > 0.0)]
    eff = hst_eff(hstdata, mcloud)
    avir = hst_alphavir(hstdata)
    tstarmax = time[np.argmax(mstar > 0.8*eff)]
    tstarplus = time[np.argmax(mstar > 0.9*eff)] - tstarmax
    tstarminus = tstarmax - time[np.argmax(mstar > 0.7*eff)]
    talphamax = time[np.argmax(avir > 5)]
    talphaplus = time[np.argmax(avir > 10)] - talphamax
    talphaminus = talphamax - time[np.argmax(avir > 2)]

    modelsigma[i] = sigmacloud
    tstars[i] = tstar
    tsmax[i] = tstarmax
    tsplus[i] = tstarplus
    tsminus[i] = tstarminus
    tamax[i] = talphamax
    taplus[i] = talphaplus
    taminus[i] = talphaminus
    tconvMyr[i] = tMyr / tff
    effs[i] = eff

plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.errorbar(modelsigma,tsmax,yerr=[tsminus,tsplus],fmt='k.')
plt.errorbar(modelsigma,tamax,yerr=[taminus,taplus],fmt='r.')
plt.text(10**(1+0.05*2),0.9*3,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,3])
plt.yticks([0,1,2,3])
#plt.xticks([1.0e1,1.0e2,1.0e3])
plt.xticks([])

plt.ylabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

axbeta = plt.subplot(2,2,2)
plt.errorbar(modelsigma,tsmax-tstars,yerr=[tsminus,tsplus],fmt='k.')
plt.errorbar(modelsigma,tamax-tstars,yerr=[taminus,taplus],fmt='r.')
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axbeta.yaxis.set_label_position("right")
plt.text(10**(1+0.05*2),0.9*2,r"$\displaystyle(b)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,2])
plt.yticks([0,1,2])
#plt.xticks([1.0e1,1.0e2,1.0e3])
plt.xticks([])

plt.ylabel(r"$\displaystyle (t - t_*) / t_{\rm ff}$")
#plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(2,2,3)
plt.errorbar(modelsigma,(tsmax-tstars)/tconvMyr,yerr=[tsminus/tconvMyr,tsplus/tconvMyr],fmt='k.')
plt.errorbar(modelsigma,(tamax-tstars)/tconvMyr,yerr=[taminus/tconvMyr,taplus/tconvMyr],fmt='r.')
plt.text(10**(1+0.05*2),0.9*10,r"$\displaystyle(c)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,10])
plt.yticks([0,5,10])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle (t - t_*) / {\rm Myr}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

axchi = plt.subplot(2,2,4)
plt.errorbar(modelsigma,(tsmax-tstars)/(tconvMyr * effs),yerr=[tsminus/(tconvMyr*effs),tsplus/(tconvMyr*effs)],fmt='k.')
plt.errorbar(modelsigma,(tamax-tstars)/(tconvMyr * effs),yerr=[taminus/(tconvMyr*effs),taplus/(tconvMyr*effs)],fmt='r.')
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axchi.yaxis.set_label_position("right")
plt.text(10**(1+0.05*2),0.9*40,r"$\displaystyle(d)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,40])
plt.yticks([0,10,20,30,40])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle t_{\rm dep} / {\rm Myr}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'f26.pdf')
pp.savefig()
pp.close()


