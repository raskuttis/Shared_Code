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
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_Ds0.3/', 'UV_M5.0e4_R15.0_N256_Tf4_Ds0.5/',
          'UV_M5.0e4_R15.0_N256_Tf4_Ds0.8/', 'UV_M5.0e4_R15.0_N256_Tf4_Ds1.0/']
nds = len(dflist)
times = []
masses = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    times.append(time)
    masses.append(mstar)

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

colors = ['k','r','b','g','c']
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
plt.subplot(1,1,1)
pa, pb, pc, pd = plt.plot(times[0],masses[0],colors[0],times[1],masses[1],colors[1],times[2],masses[2],colors[2],times[3],masses[3],colors[3])
plt.legend((pa, pb, pc, pd), (r"$\displaystyle \delta_{\rm sol} = 3/10$",r"$\displaystyle \delta_{\rm sol} = 1/2$",r"$\displaystyle \delta_{\rm sol} = 4/5$",r"$\displaystyle \delta_{\rm sol} = 1$"),prop={'size':8})
plt.axis([0,3,0,1])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'fr3e.pdf')
pp.savefig()
pp.close()


