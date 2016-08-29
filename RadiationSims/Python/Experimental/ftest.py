from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_math import *

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

datafolder = 'UV_M5.0e3_R5.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
msol = out_msol(outlines)

time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
pfit, bfit = mstar_plfit(time, mstar, 0.01, 0.3)
tcut, afpl, aspl, anfpl, anspl = mstar_brokenplfit(time, mstar, 0.01, 0.3)
tcutlin, anflin, anslin = mstar_brokenlinfit(time, mstar, 0.01, 0.3)
tstar = time[np.argmax(mstar > 0.0)]
mstarfit = broken_pl(time, tstar, 4.0, pfit, 1.0, 10**bfit, 1.0)
#mstarfit = 10**(pfit * np.log10(time-tstar) + bfit)
plmstarfit = broken_pl(time, tstar, tcut, afpl, aspl, anfpl, anspl)
linmstarfit = broken_pl(time, tstar, tcutlin, 1.0, 1.0, anflin, anslin)

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.plot(time, mstar, 'k',time,mstarfit,'r--')
plt.axis([0,2,0,0.5])
plt.xticks([])
plt.yticks([0,0.5])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(a)$")

plt.subplot(2,2,2)
plt.plot(time, mstar, 'k', time, plmstarfit, 'r--')
plt.axis([0,2,0,0.5])
plt.xticks([])
plt.yticks([])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(b)$")

plt.subplot(2,2,3)
plt.plot(time, mstar, 'k', time, linmstarfit, 'r--')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2])
plt.yticks([0,0.5])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(c)$")

plt.subplot(2,2,4)
pmstar, = plt.plot(time, mstar, 'k')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2])
plt.yticks([])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(d)$")

pp = PdfPages(plotdir + 'f23test.pdf')
pp.savefig()
pp.close()


