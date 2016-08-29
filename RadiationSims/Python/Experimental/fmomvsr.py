from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_fluxes import *
from ..Hyperion.hyp_math import *
from ..Hyperion.hyp_star import *
import matplotlib.pyplot as plt
import scipy.integrate as spint

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hostname = 'raskutti@bellona.astro.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hstfile = 'id0/RadParGrav.hst'

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_Fluxes/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
tff = out_tff(outlines)
sigma = out_vturb(outlines)
psi = out_Psi(outlines)
clight = out_clight(outlines)

pnorm = mcloud * sigma
mnorm = mcloud
vnorm = sigma

effstar = hst_mstar(hstdata, mcloud)
time = hst_time(hstdata, tff)
pstar = effstar * mcloud * psi / clight
pstar = spint.cumtrapz(pstar, time * tff, initial=0)
ke = hst_ke(hstdata)
mgas = hst_mgas(hstdata, mcloud)
pturb = np.sqrt(2.0 * mgas * mcloud * ke)

timef, frint = fluxintvst(fluxdata, 6, 1.0, rcloud, tff)
prad = 4.0 * np.pi * rcloud**2 * frint / clight
timef, mflux = fluxintvst(fluxdata, 4, 1.0, rcloud, tff)
timef, pr = fluxintvst(fluxdata, 11, 1.0, rcloud, tff)
timef, prej = fluxintvst(fluxdata, 12, 1.0, rcloud, tff)

timef, frint = fluxintvst(fluxdata, 6, 2.0, rcloud, tff)
pradout = 4.0 * np.pi * (2.0*rcloud)**2 * frint / clight
timefout, mfluxout = fluxintvst(fluxdata, 4, 2.0, rcloud, tff)
timefout, prout = fluxintvst(fluxdata, 11, 2.0, rcloud, tff)
timefout, prejout = fluxintvst(fluxdata, 12, 2.0, rcloud, tff)

timef, dprejdt = fluxvst(fluxdata, 12, 1.0, rcloud, tff)
timef, dmdt = fluxvst(fluxdata, 4, 1.0, rcloud, tff)
vr = dprejdt / dmdt

timefout, dprejdtout = fluxvst(fluxdata, 12, 2.0, rcloud, tff)
timefout, dmdtout = fluxvst(fluxdata, 4, 2.0, rcloud, tff)
vrout = dprejdtout / dmdtout

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [xsize,ysize])
plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,1,1)
pr, ps, prd, pt = plt.plot(timef,prejout/pnorm,'k',time,pstar/pnorm,'b',timef,pradout/pnorm,'r',time,pturb/pnorm,'g')
pr, ps, prd, pt = plt.plot(timef,prej/pnorm,'--k',time,pstar/pnorm,'--b',timef,prad/pnorm,'--r',time,pturb/pnorm,'--g')
plt.legend((pr, prd, ps, pt), (r"$\displaystyle p_{\rm r, ej}$",r"$\displaystyle p_{\rm *, esc}$",r"$\displaystyle p_{\rm *, in}$",r"$\displaystyle p_{\rm turb}$"),prop={'size':8},loc=2)
plt.axis([0,3,0,5.0])
plt.xticks([0,1,2,3])
plt.yticks([0,1,2,3,4,5])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle p_{\rm r} / M_{\rm cloud} \sigma$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'fmomvsr.pdf')
pp.savefig()
pp.close()




