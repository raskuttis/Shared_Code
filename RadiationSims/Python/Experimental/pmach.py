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
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4_NF/']
nds = len(dflist)
times = []
machs = []
alphas = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    csound = out_csound(outlines)
    time = hst_time(hstdata, tff)
    mass = hst_mgas(hstdata, mcloud)
    alpha = hst_alphavir(hstdata)
    machn = hst_mach(hstdata, csound)
    times.append(time)
    machs.append(machn)
    alphas.append(alpha)

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

ycap = 10**(1+0.9)
plt.subplot(1,2,1)
for i in xrange(0,nds):
    plt.plot(times[i], machs[i])
plt.axis([0,3,1.0e1,1.0e2])
plt.yscale('log')
plt.xticks([])
plt.yticks([1.0e1, 1.0e2])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \mathcal{M}$")
plt.text(0.05*3,ycap,r"$\displaystyle(a)$")

plt.subplot(1,2,2)
for i in xrange(0,nds):
    plt.plot(alphas[i], machs[i])
plt.axis([1.0e-1,1.0e2,1.0e1,1.0e2])
plt.yscale('log')
plt.xscale('log')
plt.xticks([1.0e-1,1.0,1.0e1,1.0e2])
plt.yticks([])

plt.xlabel(r"$\displaystyle \alpha_{\rm vir}$")
#plt.ylabel(r"$\displaystyle \mathcal{M}$")
plt.text(0.05*3,ycap,r"$\displaystyle(b)$")

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()


