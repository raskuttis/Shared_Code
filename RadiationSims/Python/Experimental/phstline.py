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

dflist = ['UV_M5.0e4_R15.0_N128_Tf4_B0.2/', 'UV_M5.0e4_R15.0_N128_Tf4_B0.02/', 'UV_M5.0e4_R15.0_N128_Tf4_B0.002/', 'UV_M5.0e4_R15.0_N128_Tf4_B0.0002/', 'UV_M5.0e4_R15.0_N128_Tf4_B0.00002/']
nds = len(dflist)
times = []
masses = []
num = 13

for i in xrange(0,nds):
    datafolder = dflist[i]
    print i, datafolder
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    time = hst_time(hstdata, 1.0)
    mstar = hst_num(hstdata, num)
    times.append(time)
    masses.append(mstar)

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

colors = ['k','r','b','g','c','y','m']
plt.subplot(1,1,1)
for i in xrange(0,nds):
    #plt.plot(times[i],np.log10(np.fabs(masses[i])),colors[i])
    plt.plot(times[i],masses[i],colors[i])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")

pp = PdfPages(plotdir + 'f{0:d}.pdf'.format(num))
pp.savefig()
pp.close()


