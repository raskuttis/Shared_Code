from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_star import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N128_Tf4_NF_B0.05_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_NF_B0.05_LD_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_NF_B5.0_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_NF_B5.0_LD_SDs/']
betas = [0.05, 0.05, 5.0, 5.0]
nds = len(dflist)
times = []
masses = []

pref = 959.4

for i in xrange(0,nds):
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    print max(time)
    mstar = hst_num(hstdata, 13) / mcloud
    mgas = hst_num(hstdata, 2) / mcloud
    #mstar = hst_mstar(hstdata, mcloud)
    mf = np.sqrt(pref * betas[i])
    times.append(time)
    masses.append(mstar + mgas)
    print mf

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

colors = ['k','r','b','g','c']
plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
for i in xrange(0,nds):
    plt.plot(times[i],masses[i],colors[i])
plt.axis([0,4,0,2])
plt.xticks([0,1,2,3,4])
plt.yticks([0,1.0,2.0])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'ptotmass.pdf')
pp.savefig()
pp.close()


