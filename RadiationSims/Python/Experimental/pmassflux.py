from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M5.0e4_R15.0_N128_Tf4_NF_B0.05_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_NF_B0.05_LD_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_NF_B5.0_LD_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_NF_B5.0_LD_SDs/']
betas = [0.05, 0.05, 5.0, 5.0]
nds = len(dflist)
times = []
massfluxes = []

pref = 959.4

for i in xrange(0,nds):
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    time = hst_time(hstdata, tff)
    cs = out_csound(outlines)
    gravc = out_G(outlines)
    rhobar = out_rhocloud(outlines)
    beta = betas[i]
    bz = np.sqrt(8.0 * np.pi * cs**2 * rhobar / beta)
    print max(time)
    bztote = hst_num(hstdata, 18)
    print bz, bztote[0] / ((4.0 * rcloud)**3)
    bztote = hst_num(hstdata, 17)
    bsq = 2.0 * (bztote / (4.0 * rcloud**3)) + bz**2
    bmean = np.sqrt(bsq)
    mfall = 2 * np.sqrt(gravc) * mcloud / (bmean * rcloud**2)
    print mfall[0]
    mf = np.sqrt(pref * beta * rhobar)
    times.append(time)
    massfluxes.append(mfall)

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
    plt.plot(times[i],massfluxes[i],colors[i])
plt.axis([0,4,1.0,1e3])
plt.yscale('log')
plt.xticks([])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'pmassflux.pdf')
pp.savefig()
pp.close()


