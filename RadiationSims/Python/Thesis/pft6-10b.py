from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B50.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B20.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B5.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B2.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B1.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.5_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.2_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.1_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_Fluxes/']
betas = [50.0, 20.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]
bzwrong = 1
nfs = len(dflist)

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
    rcloud = out_rcloud(outlines)
    cs = out_csound(outlines)
    gravc = out_G(outlines)
    rhobar = out_rhocloud(outlines)
    beta = betas[i]
    bz = np.sqrt(8.0 * np.pi * cs**2 * rhobar / beta)
    if bzwrong == 1:
        bz = bz * np.sqrt(4.0 * np.pi)
    mf0 = 2 * np.sqrt(gravc) * mcloud / (bz * rcloud**2)
    
    time = hst_time(hstdata, tff)
    mstar = hst_num(hstdata, 13) / mcloud
    tstar = time[np.argmax(mstar > 0.0)]
    eff = np.max(mstar)
    ke = hst_ke(hstdata)
    pe = hst_num(hstdata, 12)
    avir = -2.0 * ke / pe
    tstarmax = time[np.argmax(mstar > 0.8*eff)]
    tstarplus = time[np.argmax(mstar > 0.9*eff)] - tstarmax
    tstarminus = tstarmax - time[np.argmax(mstar > 0.7*eff)]
    talphamax = time[np.argmax(avir > 5)]
    talphaplus = time[np.argmax(avir > 10)] - talphamax
    talphaminus = talphamax - time[np.argmax(avir > 2)]

    modelsigma[i] = mf0
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
plt.axis([0.5,1.0e2,0,3])
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
plt.axis([0.5,1.0e2,0,2])
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
plt.axis([0.5,1.0e2,0,10])
plt.yticks([0,5,10])
plt.xticks([1.0,1.0e1,1.0e2])

plt.ylabel(r"$\displaystyle (t - t_*) / {\rm Myr}$")
plt.xlabel(r"$\displaystyle \mu_{\Phi,0}$")

axchi = plt.subplot(2,2,4)
plt.errorbar(modelsigma,(tsmax-tstars)/(tconvMyr * effs),yerr=[tsminus/(tconvMyr*effs),tsplus/(tconvMyr*effs)],fmt='k.')
plt.errorbar(modelsigma,(tamax-tstars)/(tconvMyr * effs),yerr=[taminus/(tconvMyr*effs),taplus/(tconvMyr*effs)],fmt='r.')
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axchi.yaxis.set_label_position("right")
plt.text(10**(1+0.05*2),0.9*40,r"$\displaystyle(d)$")
plt.xscale('log')
plt.axis([0.5,1.0e2,0,40])
plt.yticks([0,10,20,30,40])
plt.xticks([1.0,1.0e1,1.0e2])

plt.ylabel(r"$\displaystyle t_{\rm dep} / {\rm Myr}$")
plt.xlabel(r"$\displaystyle \mu_{\Phi,0}$")

pp = PdfPages(plotdir + 'ft6-10b.pdf')
pp.savefig()
pp.close()


