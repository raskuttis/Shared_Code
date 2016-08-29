from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_fluxes import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hostname = 'raskutti@bellona.astro.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r0.dat'

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_Fluxes/'
#outlines = read_outfile(hostname,datadir + datafolder + outfile)
fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)

#mcloud = out_mcloud(outlines)
#rcloud = out_rcloud(outlines)
#tff = out_tff(outlines)
tff = 4.4
rcloud = 15.0
mcloud = 1.445e6
fourpiG = 1.870283e-03
gravc = fourpiG / (4.0 * np.pi)
rad, tau = fluxvsr(fluxdata, 3, 0.0, rcloud, tff)
rad, rho = fluxvsr(fluxdata, 5, 0.0, rcloud, tff)
rad, mass = fluxintvsr(fluxdata, 5, 0.0, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, 0.0, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, 0.0, rcloud, tff)
rad, phi = fluxvsr(fluxdata, 13, 1.5, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, 1.5, rcloud, tff)
chi = rhofr / (rho * fr)
phinorm = gravc * mcloud / rcloud
dphidrnorm = phinorm / rcloud
phith = phinorm / 2.0 * (rad**2 - 3.0)
dphidrth = dphidrnorm * rad

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,2,1)
plt.plot(rad,phi/phinorm,'k',rad,phith/phinorm,'r')
plt.axis([0,3,-2.0,0.0])
plt.xticks([0,1,2,3])
plt.yticks([-2.0,-1.5,-1,-0.5,0.0])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle \Phi / (G M_{\rm cloud} / r_{\rm cloud})$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

plt.subplot(1,2,2)
plt.plot(rad,dphidr/dphidrnorm,'k',rad,dphidrth/dphidrnorm,'r')
plt.axis([0,3,0.0,1.0])
plt.xticks([0,1,2,3])
plt.yticks([0.0,0.5,1.0])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle d\Phi / dr / (G M_{\rm cloud} / r^2_{\rm cloud})$")

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()




