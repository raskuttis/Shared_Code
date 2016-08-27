from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_out import *
from hyp_fluxes import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hostname = 'raskutti@bellona.astro.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r0.dat'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_Fluxes/'
#outlines = read_outfile(hostname,datadir + datafolder + outfile)
fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)

#mcloud = out_mcloud(outlines)
#rcloud = out_rcloud(outlines)
#tff = out_tff(outlines)
tff = 4.4
rcloud = 15.0
mcloud = 1.445e6
rad, tau = fluxvsr(fluxdata, 3, 0.0, rcloud, tff)
rad, rho = fluxvsr(fluxdata, 5, 0.0, rcloud, tff)
rad, mass = fluxintvsr(fluxdata, 5, 0.0, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, 0.0, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, 0.0, rcloud, tff)
rad, phi = fluxvsr(fluxdata, 13, 0.0, rcloud, tff)
chi = rhofr / (rho * fr)

print rad
print phi

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,1,1)
plt.plot(rad,mass/mcloud,'k')
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1,1.5])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()




