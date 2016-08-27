from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_out import *
from hyp_hst import *
from hyp_fluxes import *
from hyp_math import *
from hyp_star import *
from hyp_pdf import *
from hyp_models import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hostname = 'raskutti@tiger.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hstfile = 'id0/RadParGrav.hst'
starfile = 'star'
nconv = 1
nconvpdf = 25

dflist = set_model_lookup('Fluxes/')
dfnflist = set_model_lookup('NF_Fluxes/')
nfs = len(dflist)

nts = 100
plotts = np.linspace(0.0,3.0,num=nts)
funbs = np.zeros(nfs)
feunbs = np.zeros(nfs)
allsigmas = np.zeros(nfs)
qval = 0.68

for i in xrange(2,nfs):
    
    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    clight = out_clight(outlines)
    tff = out_tff(outlines)
    psi = out_Psi(outlines)
    kappa = out_kappa(outlines)
    gravc = out_G(outlines)
    rhocloud = out_rhocloud(outlines)
    sigmacloud = out_sigma(outlines)
    allsigmas[i] = sigmacloud

    time = hst_time(hstdata, tff)
    eff = hst_mstar(hstdata, mcloud)
    mgas = hst_mgas(hstdata, mcloud)
    eps = hst_eff(hstdata, mcloud)
    stesc = np.max(eff) - eff[-1]
    stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)
    funb = star_funboundattt(stardata, time, 4.0, gravc, mcloud)
    funbs[i] = (funb + stesc) / eps
    feunbs[i] = funb + stesc
    print i, funbs[i], feunbs[i], funb, stesc, eps

fnfunbs = np.zeros(nfs)
fenfunbs = np.zeros(nfs)
qval = 0.68

for i in xrange(2,nfs):
    
    datafolder = dfnflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    clight = out_clight(outlines)
    tff = out_tff(outlines)
    psi = out_Psi(outlines)
    kappa = out_kappa(outlines)
    gravc = out_G(outlines)
    rhocloud = out_rhocloud(outlines)
    sigmacloud = out_sigma(outlines)
    
    time = hst_time(hstdata, tff)
    eff = hst_mstar(hstdata, mcloud)
    mgas = hst_mgas(hstdata, mcloud)
    eps = hst_eff(hstdata, mcloud)
    stesc = np.max(eff) - eff[-1]
    stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)
    funb = star_funboundattt(stardata, time, 4.0, gravc, mcloud)
    fnfunbs[i] = (funb + stesc) / eps
    fenfunbs[i] = funb + stesc
    print i, fnfunbs[i], fenfunbs[i], funb, stesc, eps

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
pa, pb = plt.plot(allsigmas,funbs,'k.',allsigmas,fnfunbs,'r.')
plt.legend((pa, pb), (r"$\displaystyle {\rm Feedback}$",r"$\displaystyle {\rm No Feedback}$"),prop={'size':8})
plt.text(10**(1+0.05*2),0.9,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,0.5])
plt.yticks([0,0.25,0.5])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle f_{\rm unb}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(1,2,2)
plt.subplots_adjust(bottom=0.2)
plt.plot(allsigmas,feunbs,'k.',allsigmas,fenfunbs,'r.')
plt.text(10**(1+0.05*2),0.9,r"$\displaystyle(b)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,0.5])
plt.yticks([0,0.25,0.5],[' ',' ',' '])
plt.xticks([1.0e1,1.0e2,1.0e3])

#plt.ylabel(r"$\displaystyle f_{\rm unb}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'fbounds.pdf')
pp.savefig()
pp.close()




