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
fluxfile = 'fluxes_r0.dat'
hstfile = 'id0/RadParGrav.hst'

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_Fluxes/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
clight = out_clight(outlines)
tff = out_tff(outlines)
psi = out_Psi(outlines)
kappa = out_kappa(outlines)
nconv = 10

tplot = 1.0
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
fedd = rhofr * kappa / (rhodphidr * clight)
feddspec = fr * kappa / (dphidr * clight)
fedd = np.convolve(fedd, np.ones((nconv,))/nconv, mode='same')
feddspec = np.convolve(feddspec, np.ones((nconv,))/nconv, mode='same')

tplot = 1.5
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
feddlate = rhofr * kappa / (rhodphidr * clight)
feddspeclate = fr * kappa / (dphidr * clight)
feddlate = np.convolve(feddlate, np.ones((nconv,))/nconv, mode='same')
feddspeclate = np.convolve(feddspeclate, np.ones((nconv,))/nconv, mode='same')

timef, rhofr = fluxvst(fluxdata, 7, 1.0, rcloud, tff)
timef, rhodphidr = fluxvst(fluxdata, 15, 1.0, rcloud, tff)
timef, fr = fluxvst(fluxdata, 6, 1.0, rcloud, tff)
timef, dphidr = fluxvst(fluxdata, 14, 1.0, rcloud, tff)
feddf = rhofr * kappa / (rhodphidr * clight)
feddfspec = fr * kappa / (dphidr * clight)

timef, rhofr = fluxvst(fluxdata, 7, 2.0, rcloud, tff)
timef, rhodphidr = fluxvst(fluxdata, 15, 2.0, rcloud, tff)
timef, fr = fluxvst(fluxdata, 6, 2.0, rcloud, tff)
timef, dphidr = fluxvst(fluxdata, 14, 2.0, rcloud, tff)
feddfout = rhofr * kappa / (rhodphidr * clight)
feddfspecout = fr * kappa / (dphidr * clight)

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
f, fs = plt.plot(rad,fedd,'k',rad,feddlate,'r')
plt.legend((f, fs), (r"$\displaystyle t = t_{\rm ff}$",r"$\displaystyle t = 1.5~t_{\rm ff}$"),prop={'size':8},loc=2)
plt.axis([0,2,0,10.0])
plt.xticks([0,1,2])
plt.yticks([0,5,10])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle f_{\rm edd}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

plt.subplot(1,2,2)
f, fs = plt.plot(timef,feddf,'k',timef,feddfout,'r')
plt.legend((f, fs), (r"$\displaystyle r = r_{\rm cloud}$",r"$\displaystyle r = 2~r_{\rm cloud}$"),prop={'size':8},loc=2)
plt.axis([0,2,0,10.0])
plt.xticks([0,1,2])
plt.yticks([0,5,10],[' ',' ',' '])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle f_{\rm abs}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'ffedd.pdf')
pp.savefig()
pp.close()




