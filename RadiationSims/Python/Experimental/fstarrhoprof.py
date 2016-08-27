from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_out import *
from hyp_hst import *
from hyp_fluxes import *
from hyp_math import *
from hyp_star import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hostname = 'raskutti@tiger.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hstfile = 'id0/RadParGrav.hst'
starfile = 'star'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_Fluxes/'
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

time = hst_time(hstdata, tff)
eff = hst_mstar(hstdata, mcloud)
mgas = hst_mgas(hstdata, mcloud)
stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)
nconv = 20

tplot = 0.57
allstarrad, allstarmvsr = star_massvsr(stardata)
tmini = np.abs(time - tplot).argmin()
starrad = allstarrad[tmini,:]
starmvsr = allstarmvsr[tmini,:]
sortinds = np.argsort(starrad)
starrad = starrad[sortinds]
starmvsr = starmvsr[sortinds]
starmvsr = np.cumsum(starmvsr)
rstarearly = starrad / rcloud
mstarearly = starmvsr / mcloud

tplot = 1.06
tmini = np.abs(time - tplot).argmin()
starrad = allstarrad[tmini,:]
starmvsr = allstarmvsr[tmini,:]
sortinds = np.argsort(starrad)
starrad = starrad[sortinds]
starmvsr = starmvsr[sortinds]
starmvsr = np.cumsum(starmvsr)
rstarout = starrad / rcloud
mstarout = starmvsr / mcloud

tplot = 1.55
tmini = np.abs(time - tplot).argmin()
starrad = allstarrad[tmini,:]
starmvsr = allstarmvsr[tmini,:]
sortinds = np.argsort(starrad)
starrad = starrad[sortinds]
starmvsr = starmvsr[sortinds]
starmvsr = np.cumsum(starmvsr)
rstarmid = starrad / rcloud
mstarmid = starmvsr / mcloud

tplot = 2.07
tmini = np.abs(time - tplot).argmin()
starrad = allstarrad[tmini,:]
starmvsr = allstarmvsr[tmini,:]
sortinds = np.argsort(starrad)
starrad = starrad[sortinds]
starmvsr = starmvsr[sortinds]
starmvsr = np.cumsum(starmvsr)
rstarlate = starrad / rcloud
mstarlate = starmvsr / mcloud

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
t1, t2, t3, t4 = plt.plot(rstarearly,mstarearly,'k',rstarout,mstarout,'r',rstarmid,mstarmid,'b',rstarlate,mstarlate,'g')
plt.legend((t1,t2,t3,t4), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of,final}$"),prop={'size':8})
plt.axis([1.0e-1,2.0,1.0e-4,1.0e2])
plt.yscale('log')
plt.yticks([1.0e-4,1.0e-2,1.0,1.0e2])
plt.xscale('log')
plt.ylabel(r"$\displaystyle M_*(r) / M_{\rm cl,0}$")
plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")

pp = PdfPages(plotdir + 'flogstarrhoprof.pdf')
pp.savefig()
pp.close()




