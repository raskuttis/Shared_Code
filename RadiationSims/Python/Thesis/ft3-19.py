from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_models import *
import matplotlib.pyplot as plt

## Plot showing evolution of stellar mass, gas mass, outflow mass and mach number for model with no gravity

## Define the locations where data is located and where to plot from ..Hyperion.hyp_models
hostname, datadir, hstfile, outfile, plotdir = init_dirs('ast')
fname = 'ft3-19.pdf'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_NG/'
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

print alphavir

## Plotting setup
ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
pmstar, pmgas, pmof, = plt.plot(time, mstar, 'k', time, mgas, 'r', time, mof, 'b', linewidth = 1.0)
plt.legend((pmstar, pmgas, pmof), (r"$\displaystyle M_*$",r"$\displaystyle M_{\rm gas}$",r"$\displaystyle M_{\rm of}$"),prop={'size':8})
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1,1.5])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

#axmach = axalpha.twinx()
axmach = plt.subplot(1,2,2)
axmach.yaxis.set_label_position("right")
axmach.plot(time, machn,'k')
plt.yscale('log')
plt.axis([0,3,5.0,50])
plt.yticks([10.0])
plt.xticks([0,1,2,3])

plt.ylabel(r"$\displaystyle \mathcal{M}$")
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle \alpha_{\rm vir}$")
plt.text(0.05*3,10**(0.9+0.7),r"$\displaystyle(b)$")

pp = PdfPages(plotdir + fname)
pp.savefig()
pp.close()


