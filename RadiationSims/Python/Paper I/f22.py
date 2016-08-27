from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_star import *
from hyp_math import *
import numpy as np
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4_NF/']
nds = len(dflist)
times = []
masses = []
massfits = []
tfits = np.logspace(-2.0, 1.0, num=100, base=10.0)

for i in xrange(0,nds):
    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    tsind = np.argmax(mstar > 0.0)
    tstar = time[tsind]
    times.append(time - tstar)
    masses.append(mstar)
    pfit, bfit = mstar_plfit(time, mstar, 0.1, 0.3)
    mstarfit = 10**(pfit * np.log10(tfits) + bfit)
    massfits.append(mstarfit)


plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,1,1)
plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
plt.plot(times[0], masses[0], 'k', times[1], masses[1], 'r',tfits,massfits[0],'k--',tfits,massfits[1],'r--')

plt.axis([0.04,4.0,1.0e-2,1.0])
plt.xscale('log')
plt.yscale('log')
plt.xticks([0.1,1.0],['0.1','1.0'])
plt.yticks([1.0e-2,1.0e-1,1])

plt.xlabel(r"$\displaystyle (t - t_*) / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'f30.pdf')
pp.savefig()
pp.close()


