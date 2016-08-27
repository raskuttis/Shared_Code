from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_math import *
from hyp_pdf import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'denpdfall.dat'
starfile = 'star'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

datafolder = 'UV_M2.0e5_R25.0_N256_Tf4_NF/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
nconv = 50

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
rhocloud = out_rhocloud(outlines)
msol = out_msol(outlines)
time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
#mstar = np.convolve(mstar, np.ones((nconv,))/nconv, mode='same')
tstar = time[np.argmax(mstar > 0.0)]
mgas = hst_mgas(hstdata, mcloud)
eff = hst_eff(hstdata, mcloud)
maxm = min([0.2,0.9*eff])

pfit, bfit = mstar_plfit(time, mstar, 0.01, maxm)
tcut, afpl, aspl, anfpl, anspl = mstar_brokenplfit(time, mstar, 0.01, maxm)
tcutlin, anflin, anslin = mstar_brokenlinfit(time, mstar, 0.01, maxm)
tstar = time[np.argmax(mstar > 0.0)]
mstarfit = broken_pl(time, tstar, 4.0, pfit, 1.0, 10**bfit, 1.0)
mstarfit = 10**(pfit * np.log10(time-tstar) + bfit)
plmstarfit = broken_pl(time, tstar, tcut, afpl, aspl, anfpl, anspl)
linmstarfit = broken_pl(time, tstar, tcutlin, 1.0, 1.0, anflin, anslin)
mcut = mstar[np.argmax(time > tcut)]

print tcutlin, tcut, mcut, aspl, anflin, anslin

# Read in PDF for calculation of mean density
pdftime = sim_pdftime(pdflines, tff)
nlines = len(pdflines)
nts = int(np.floor(nlines / 2))
times = []
besteff = 0.45
besteffbreak = 1.21

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.plot(time, mstar, 'k',time,mstarfit,'r--')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2],[' ',' ',' '])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(a)$")

plt.subplot(2,2,2)
plt.plot(time, mstar, 'k', time, plmstarfit, 'r--')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2],[' ',' ',' '])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5],[' ',' ',' ',' ',' ',' '])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(b)$")

plt.subplot(2,2,3)
plt.plot(time, mstar, 'k', time, linmstarfit, 'r--')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(c)$")

plt.subplot(2,2,4)
plt.plot(time, mstar, 'k', time, linmstarfit, 'r--')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5],[' ',' ',' ',' ',' ',' '])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(d)$")

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()


