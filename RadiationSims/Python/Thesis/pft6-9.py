from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_star import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B5.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_SDs/']
dftlist = [1,0,0]
nds = len(dflist)
nconv = 100
times = []
tsfrs = []
masses = []
sfrs = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    dftype = dftlist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    if dftype == 1:
        mstar = hst_mstar(hstdata, mcloud)
        mgas = hst_mgas(hstdata, mcloud)
        eff = hst_eff(hstdata, mcloud)
    else:
        mstar = hst_num(hstdata, 13) / mcloud
        mgas = hst_num(hstdata, 2) / mcloud
        eff = np.max(mstar)
    msconv = np.convolve(mstar, np.ones((nconv,))/nconv, mode='same')
    sfr = np.diff(msconv) / np.diff(time)
    tsfr = time[1:]
    tsfrs.append(tsfr)
    times.append(time)
    masses.append(mstar)
    sfrs.append(np.convolve(sfr, np.ones((nconv,))/nconv, mode='same'))

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
pa, pb, pc, = plt.plot(times[0], masses[0], 'k', times[1], masses[1], 'r', times[2], masses[2], 'b')
plt.legend((pa, pb, pc), (r"$\displaystyle {\rm Run~I}$",r"$\displaystyle {\rm Run~MIII}$",r"$\displaystyle {\rm Run~MI}$"),prop={'size':8},loc=1)
plt.axis([0,3,0,1])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \varepsilon$")
plt.text(0.9*3,0.9*1,r"$\displaystyle(a)$")

axsfr = plt.subplot(1,2,2)
axsfr.yaxis.set_label_position("right")
plow, pmid, phigh, = plt.plot(tsfrs[0], sfrs[0], 'k', tsfrs[1], sfrs[1], 'r', tsfrs[2], sfrs[2], 'b')
plt.axis([0,3,0.0,1.0])
#plt.yscale('log')
plt.xticks([0,1,2,3])
#plt.yticks([1.0e-2,1.0e-1,1.0])
plt.yticks([0.0,0.5,1.0])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \varepsilon_{\rm ff}$")
#plt.text(0.05*3,10**(-1+0.9*3),r"$\displaystyle(b)$")
plt.text(0.05*3,0.9,r"$\displaystyle(b)$")

pp = PdfPages(plotdir + 'ft6-9.pdf')
pp.savefig()
pp.close()


