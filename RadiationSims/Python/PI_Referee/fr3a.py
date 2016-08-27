from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_star import *

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

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,1,1)
plt.subplots_adjust(left=0.2)
plt.subplots_adjust(right=0.8)
pmstar, pmgas, pmof, = plt.plot(time, mstar, 'k', time, mgas, 'r', time, mof, 'b', linewidth = 1.0)
plt.legend((pmstar, pmgas, pmof), (r"$\displaystyle M_*$",r"$\displaystyle M_{\rm gas}$",r"$\displaystyle M_{\rm of}$"))
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3],[' ',' ',' ',' '])
plt.yticks([0,0.5,1,1.5])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

#axalpha = plt.subplot(2,1,2)
#plt.subplots_adjust(wspace=0.3)
#palphavir = plt.plot(time, alphavir, 'k', linewidth = 1.0)
#plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
#axalpha.yaxis.set_label_position("right")
#plt.yscale('log')
#plt.axis([0,3,0.1,100])
#plt.xticks([0,1,2,3])
#plt.yticks([0.1,1,10,100])

#axmach = axalpha.twinx()
axmach = plt.subplot(2,1,2)
axmach.plot(time, machn,'k')
plt.yscale('log')
plt.axis([0,3,5.0,50])
plt.yticks([10.0])
plt.xticks([0,1,2,3])

plt.ylabel(r"$\displaystyle \mathcal{M}$")
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle \alpha_{\rm vir}$")
plt.text(0.05*3,10**(0.9+0.7),r"$\displaystyle(b)$")

pp = PdfPages(plotdir + 'fr3b.pdf')
pp.savefig()
pp.close()


