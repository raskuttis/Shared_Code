from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_out import *
from hyp_hst import *
from hyp_fluxes import *
from hyp_math import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hostname = 'raskutti@bellona.astro.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hstfile = 'id0/RadParGrav.hst'

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_Fluxes/'
print 'Reading Out'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
print 'Reading Flux'
fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
print 'Reading Hst'
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
tff = out_tff(outlines)
psi = out_Psi(outlines)
kappa = out_kappa(outlines)

tplot = 1.06
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
chi = rhofr / (rho * fr)

tplot = 1.57
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
chilate = rhofr / (rho * fr)

timef, rho = fluxvst(fluxdata, 5, 1.0, rcloud, tff)
timef, fr = fluxvst(fluxdata, 6, 1.0, rcloud, tff)
timef, rhofr = fluxvst(fluxdata, 7, 1.0, rcloud, tff)
chif = rhofr / (rho * fr)
timef, rho = fluxvst(fluxdata, 5, 2.0, rcloud, tff)
timef, fr = fluxvst(fluxdata, 6, 2.0, rcloud, tff)
timef, rhofr = fluxvst(fluxdata, 7, 2.0, rcloud, tff)
chifout = rhofr / (rho * fr)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,ysize])
#plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,2,1)
pe, pl = plt.plot(rad,chi,'k',rad,chilate,'r')
plt.legend((pe, pl), (r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$"),prop={'size':8})
plt.axis([0,2,0,1.0])
plt.xticks([0,1,2])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle \chi_{\rho, F_r}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

plt.subplot(1,2,2)
pi, po = plt.plot(timef,chif,'k',timef,chifout,'r')
plt.legend((pi, po), (r"$\displaystyle r = r_{\rm cloud}$",r"$\displaystyle r = r_{\rm box}$"),prop={'size':8})
plt.axis([0,2,0,1.0])
plt.xticks([0,1,2])
plt.yticks([0,0.5,1],[' ',' ',' '])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle f_{\rm abs}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'fchirho.pdf')
pp.savefig()
pp.close()




