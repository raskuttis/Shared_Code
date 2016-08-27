from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_star import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
denflfile = 'denfloorall.dat'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B50.0_CR100/', 'UV_M5.0e4_R15.0_N256_Tf4_B50.0_CR250/',
          'UV_M5.0e4_R15.0_N256_Tf4_B50.0_CR500/', 'UV_M5.0e4_R15.0_N256_Tf4_B50.0_CR1000/',
          'UV_M5.0e4_R15.0_N256_Tf4_B0.05_CR100/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_CR250/',
          'UV_M5.0e4_R15.0_N256_Tf4_B0.05_CR500/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_CR1000/']
nds = len(dflist)
times = []
masses = []
dfltimes = []
dflmasses = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

    dx = out_dx(outlines)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_num(hstdata, 13) / mcloud
    
    times.append(time)
    masses.append(mstar)

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

colors = ['k','r','b','g','c']
plt.subplots_adjust(bottom=0.2)
plt.subplot(1,2,1)
pa, pb, pc, pd = plt.plot(times[0],masses[0],colors[0],times[1],masses[1],colors[1],times[2],masses[2],colors[2],times[3],masses[3],colors[3])
plt.axis([0,3,0,1])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1,r"$\displaystyle(a): \mu_{\Phi,0} \rightarrow \infty$")

axd = plt.subplot(1,2,2)
axd.yaxis.set_label_position("right")
pa, pb, pc, pd = plt.plot(times[4],masses[4],colors[0],times[4],masses[4],colors[1],times[4],masses[4],colors[2],times[4],masses[4],colors[3])
plt.legend((pa, pb, pc, pd), (r"$\displaystyle \hat{c} = 100~{\rm km~s^{-1}}$",r"$\displaystyle 250~{\rm km~s^{-1}}$",r"$\displaystyle 500~{\rm km~s^{-1}}$",r"$\displaystyle 1000~{\rm km~s^{-1}}$"),prop={'size':8})
plt.axis([0,3,0.1,1.0])
plt.xticks([0,1,2,3])
plt.yticks([0.0,0.5,1.0],[' ',' ',' '])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle M_{\rm add} / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1,r"$\displaystyle(b): \mu_{\Phi,0} = 3$")

pp = PdfPages(plotdir + 'ft3-26.pdf')
pp.savefig()
pp.close()


