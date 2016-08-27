from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_fluxes import *
from hyp_math import *
from hyp_models import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.special import erf
from scipy.stats import lognorm

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r0.dat'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B50.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B20.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B5.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B2.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B1.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.5_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.2_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.1_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_Fluxes/']
betas = [50.0, 20.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]
bzwrong = 1
nfs = len(dflist)

modelsigma = []
modeleff = []
modeleffof = []
psicgs = 2000.0
sigmaedd = surf_edd(psicgs)

altfabs = np.zeros(nfs)
allfabs = np.zeros(nfs)
allsigma = np.zeros(nfs)
allfedd = np.zeros((3,nfs))
allvinf = np.zeros((nfs, 3))
altvinf = np.zeros((nfs, 3))
alleps = np.zeros(nfs)

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    psi = out_Psi(outlines)
    clight = out_clight(outlines)
    kappa = out_kappa(outlines)
    rcloud = out_rcloud(outlines)
    vturb = out_vturb(outlines)
    vesc = vturb * np.sqrt(5.0 * 2.0 / (3.0 * 2.0))
    gravc = out_G(outlines)
    eps = hst_eff(hstdata, mcloud)
    msol = out_msol(outlines)
    cs = out_csound(outlines)
    gravc = out_G(outlines)
    rhobar = out_rhocloud(outlines)
    beta = betas[i]
    bz = np.sqrt(8.0 * np.pi * cs**2 * rhobar / beta)
    if bzwrong == 1:
        bz = bz * np.sqrt(4.0 * np.pi)
    mf0 = 2 * np.sqrt(gravc) * mcloud / (bz * rcloud**2)
    
    time = hst_time(hstdata, tff)
    eff = hst_num(hstdata, 13) / mcloud
    eps = np.max(eff)
    alleps[i] = eps
    #pnorm = mcloud * vturb
    pnorm = eps * mcloud * vesc
    pnormalt = mcloud * (1.0 - eps) * vesc
    print i, datafolder, max(time), sigmacloud, pnorm
    fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
    
    timefout, prejout = fluxintvst(fluxdata, 12, 2.0, rcloud, tff)
    
    allsigma[i] = mf0
    allfabs[i] = np.max(prejout) / pnorm
    altfabs[i] = np.max(prejout) / pnormalt

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

axv = plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
psim = plt.plot(allsigma,altfabs,'k.')
plt.xscale('log')
plt.axis([0.5,1.0e2,0,3.0])
plt.yticks([0,1,2,3])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle p_{\rm r} / (M_{\rm of} v_{\rm esc})$")
plt.xlabel(r"$\displaystyle \mu_{\Phi,0}$")

pp = PdfPages(plotdir + 'ft6-10c.pdf')
pp.savefig()
pp.close()

