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

dflist = ['UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4_NF/',
          'UV_M5.0e4_R15.0_N256_Tf4_NG_RST/', 'UV_M5.0e4_R15.0_N256_Tf4_NG_NF_RST/']
nds = len(dflist)
times = []
masses = []
mofs = []
effs = []

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
    eff = hst_eff(hstdata, mcloud)
    times.append(time)
    masses.append(mstar)
    mofs.append(mof)
    effs.append(eff)

print effs
print effs[0] / (effs[1])
print effs[2] / (effs[3])

plt.figure(figsize = [xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(2,1,1)
plt.subplots_adjust(left=0.2)
colors = ['k','r','k--','r--']
for i in xrange(0,nds):
    plt.plot(times[i], masses[i],colors[i])
plt.axis([0,3,0,1])
plt.xticks([0,1,2,3],[' ',' ',' ',' '])
plt.yticks([0,0.5,1])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M_* / M_{\rm cl,0}$")
plt.text(0.05*3,0.9,r"$\displaystyle(a)$")

plt.subplot(2,1,2)
for i in xrange(0,nds):
    plt.plot(times[i], mofs[i],colors[i])
plt.axis([0,3,0,1])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M_{\rm of} / M_{\rm cl,0}$")
plt.text(0.05*3,0.9,r"$\displaystyle(b)$")

pp = PdfPages(plotdir + 'f32.pdf')
pp.savefig()
pp.close()


