from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/',
          'UV_M5.0e4_R15.0_N256_Tf4_Cs0.5_RST/',
          'UV_M5.0e4_R15.0_N256_Tf4_Cs1.0_RST/']
nds = len(dflist)
times = []
masses = []
mofs = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    mof = hst_mofsimple(hstdata, mcloud)
    times.append(time)
    masses.append(mstar)
    mofs.append(mof)

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='10')

colors = ['k','r','b','g','c']
plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plow, pmid, phigh, = plt.plot(times[0],masses[0],colors[0],times[1],masses[1],colors[1],times[2],masses[2],colors[2])
plt.legend((plow, pmid, phigh), (r"$\displaystyle c_s = 0.2~{\rm km~s^{-1}}$",r"$\displaystyle 0.5$",r"$\displaystyle 1.0$"),prop={'size':8})
plt.axis([0,3,0,1])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M_* / M_{\rm cl,0}$")
plt.text(0.05*3,0.9,r"$\displaystyle(a)$")

axy = plt.subplot(1,2,2)
axy.yaxis.set_label_position("right")
for i in xrange(0,nds):
    plt.plot(times[i], mofs[i],colors[i])
plt.axis([0,3,0,1])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M_{\rm of} / M_{\rm cl,0}$")
plt.text(0.05*3,0.9,r"$\displaystyle(b)$")

pp = PdfPages(plotdir + 'ft3-23.pdf')
pp.savefig()
pp.close()


