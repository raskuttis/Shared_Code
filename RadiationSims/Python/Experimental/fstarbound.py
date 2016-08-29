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
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hostname = 'raskutti@bellona.astro.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hstfile = 'id0/RadParGrav.hst'
starfile = 'star'
nconv = 1
nconvpdf = 25

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
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
stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)
nts = 100
plotts = np.linspace(0,3.0,num=nts)
funbs = np.zeros(nts)
ffracs = np.zeros(nts)
for i in xrange(0,nts):
    funbs[i] = star_funboundattt(stardata, time, plotts[i], gravc, mcloud)
    tmini = np.abs(time - plotts[i]).argmin()
    effi = eff[tmini]
    ffracs[i] = funbs[i] / effi
    print i, funbs[i], ffracs[i]

rcom = star_mu(stardata) / rcloud
dcom = star_sigma(stardata) / rcloud
maxcom = star_max(stardata) / rcloud
vcom = star_vcom(stardata)

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
rc,dc,mc = plt.plot(time,rcom,'k',time,dcom,'r',time,maxcom,'b')
plt.legend((rc, dc, mc), (r"$\displaystyle \mu_*$",r"$\displaystyle \sigma_*$",r"$\displaystyle r_{max,*}$"),prop={'size':8},loc=2)
plt.axis([0,3,0,2])
plt.xticks([0,1,2,3])
plt.yticks([0,1,2])
plt.text(0.9*3,0.9*2,r"$\displaystyle(a)$")
plt.ylabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")

plt.subplot(1,2,2)
plt.subplots_adjust(bottom=0.2)
vc = plt.plot(plotts,funbs,'k',plotts,ffracs,'r')
plt.axis([0,3,0,1])
plt.xticks([0,1,2,3])
plt.text(0.9*3,0.9*1,r"$\displaystyle(b)$")
plt.ylabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")

pp = PdfPages(plotdir + 'fstarbound.pdf')
pp.savefig()
pp.close()




