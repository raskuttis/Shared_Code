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

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs0.3/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs0.4/',
          'UV_M5.0e4_R15.0_N256_Tf4_Cs0.5/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs1.5/', 'UV_M5.0e4_R15.0_N256_Tf4_Cs2.0/']
nds = len(dflist)
times = []
masses = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    mass = hst_mgas(hstdata, mcloud)
    times.append(time)
    masses.append(mass)

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,1,1)
for i in xrange(0,nds):
    plt.plot(times[i], masses[i])
plt.axis([0,3,0,1])
plt.xticks([])
plt.yticks([0,0.5,1])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()


