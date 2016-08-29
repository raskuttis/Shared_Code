from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'sdpdfxy.dat'
starfile = 'star'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

csound = out_csound(outlines)
tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
msol = out_msol(outlines)

time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
mgas = hst_mgas(hstdata, mcloud)
mof = mgas[0] - mstar - mgas
alphavir = hst_alphavir(hstdata)
machn = hst_mach(hstdata, csound)
rs = hst_rgauss(hstdata)
rshort = rs[:,2] / rcloud
rmed = rs[:,1] / rcloud
rlong = rs[:,0] / rcloud

stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)
smass = star_mass(stardata) / mcloud
smu = star_mu(stardata) / rcloud
ssigma = star_sigma(stardata) / rcloud

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
pmstar, pmgas, pmof, = plt.plot(time, mstar, 'k', time, mgas, 'r', time, mof, 'b', linewidth = 1.0)
plt.legend((pmstar, pmgas, pmof), (r"$\displaystyle M_*$",r"$\displaystyle M_{\rm gas}$",r"$\displaystyle M_{\rm of}$"))
plt.axis([0,3,0,1.5])
plt.xticks([])
plt.yticks([0,0.5,1,1.5])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

axalpha = plt.subplot(2,2,2)
plt.subplots_adjust(wspace=0.3)
palphavir = plt.plot(time, alphavir, 'k', linewidth = 1.0)
#plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
#axalpha.yaxis.set_label_position("right")
plt.yscale('log')
plt.axis([0,3,0.1,100])
plt.xticks([])
plt.yticks([0.1,1,10,100])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \alpha_{\rm vir}$")
plt.text(0.05*3,50.0,r"$\displaystyle(b)$")

axmach = axalpha.twinx()
axmach.plot(time, machn,'r')
plt.yscale('log')
plt.axis([0,3,10.0,100])
plt.yticks([10,100])

plt.ylabel(r"$\displaystyle \mathcal{M}$")

plt.subplot(2,2,3)
prshort, prmed, prlong, = plt.plot(time, rshort, 'k', time, rmed, 'r', time, rlong, 'b', linewidth = 1.0)
plt.legend((pmstar, pmgas, pmof), (r"$\displaystyle R_1$",r"$\displaystyle R_2$",r"$\displaystyle R_3$"))
plt.axis([0,3,0,4])
plt.xticks([0,1,2,3])
plt.yticks([0,1,2,3,4])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle R / R_{\rm 0}$")
plt.text(0.05*3,0.9*4,r"$\displaystyle(c)$")

axstar = plt.subplot(2,2,4)
pmu, psigma, = plt.plot(time, smu, 'k', time, ssigma, 'r', linewidth = 1.0)
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axstar.yaxis.set_label_position("right")
plt.legend((pmstar, pmgas, pmof), (r"$\displaystyle \mu_\star$",r"$\displaystyle \sigma_\star$"))
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1,1.5])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle R_\star / R_{\rm 0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(d)$")

pp = PdfPages(plotdir + 'f5.pdf')
pp.savefig()
pp.close()


