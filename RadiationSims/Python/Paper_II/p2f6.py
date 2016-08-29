from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_fluxes import *
from ..Hyperion.hyp_math import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
rhostname = 'raskutti@tiger.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hstfile = 'id0/RadParGrav.hst'
starfile = 'star'
nconv = 1
nconvpdf = 25

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M2.0e4_R15.0_N256_Tf4/', 'UV_M2.0e5_R15.0_N256_Tf4/']
dflist = ['UV_M5.0e4_R15.0_N256_Tf4_Fluxes/', 'UV_M2.0e4_R15.0_N256_Tf4_Fluxes/', 'UV_M2.0e5_R15.0_N128_Tf4_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_Fluxes/', 'UV_M2.0e4_R15.0_N256_Tf4_NF_Fluxes/', 'UV_M2.0e5_R15.0_N256_Tf4_NF_Fluxes/']
nfs = len(dflist)
nts = 100
plotts = np.linspace(0.0,3.0,num=nts)
allfs = np.zeros((nfs,nts))
allds = np.zeros((nfs,nts))
qval = 0.68

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    outlines = read_outfile(rhostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(rhostname,datadir + datafolder + hstfile)
    fluxdata = read_fluxfile(rhostname,datadir + datafolder + fluxfile)

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
    stardata = read_allstars(rhostname,datadir + datafolder + starfile,time,tff)
    print i, datafolder, max(time)

    for j in xrange(0,nts):
        allfs[i,j] = star_quartileatt(stardata, qval, time, plotts[j]) / rcloud
        rad, mr = fluxintvsr(fluxdata, 5, plotts[j], rcloud, tff)
        mr = mr / np.max(mr)
        goodind = np.argmin(np.fabs(mr - qval))
        allds[i,j] = rad[goodind]
        print j, allfs[i,j], allds[i,j]

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
rc,dc,mc = plt.plot(plotts,np.squeeze(allfs[0,:]),'k',plotts,np.squeeze(allfs[1,:]),'r',plotts,np.squeeze(allfs[2,:]),'b')
plt.plot(plotts,np.squeeze(allfs[3,:]),':k',plotts,np.squeeze(allfs[4,:]),':r',plotts,np.squeeze(allfs[5,:]),':b')
plt.legend((rc, dc, mc), (r"$\displaystyle \rm Run~I$",r"$\displaystyle \rm Run~II$",r"$\displaystyle \rm Run~III$"),prop={'size':8},loc=2)
plt.axis([0,3,0,2])
plt.xticks([0,1,2,3])
plt.yticks([0,1,2])
plt.text(0.9*3,0.9*2,r"$\displaystyle(a)$")
plt.ylabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")

plt.subplot(1,2,2)
plt.subplots_adjust(bottom=0.2)
rc,dc,mc = plt.plot(plotts,np.squeeze(allds[0,:]),'k',plotts,np.squeeze(allds[1,:]),'r',plotts,np.squeeze(allds[2,:]),'b')
plt.plot(plotts,np.squeeze(allds[3,:]),':k',plotts,np.squeeze(allds[4,:]),':r',plotts,np.squeeze(allds[5,:]),':b')
plt.axis([0,3,0,2])
plt.xticks([0,1,2,3])
plt.yticks([0,1,2],[' ',' ',' '])
plt.text(0.9*3,0.9*2,r"$\displaystyle(b)$")
#plt.ylabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")

pp = PdfPages(plotdir + 'f6.pdf')
pp.savefig()
pp.close()




