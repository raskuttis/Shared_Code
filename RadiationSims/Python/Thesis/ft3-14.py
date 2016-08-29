from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_models import *
import matplotlib.pyplot as plt

## Plot showing evolution of stellar mass and virial parameter as a function of resolution

## Define the locations where data is located and where to plot from ..Hyperion.hyp_models
hostname, datadir, hstfile, outfile, plotdir = init_dirs('ast')
fname = 'ft3-14.pdf'

## Fiducial models at different resolution (Change the last)
dflist = ['UV_M5.0e4_R15.0_N128_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4/', 'UV_M5.0e4_R15.0_N256_Tf4_L1.5/']
nds = len(dflist)
times = []
masses = []
alphas = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    ## This is a hack - Change this
    if i == nds - 1:
        oldeps = max(mstar)
    mstar = hst_mstar(hstdata, mcloud)
    if i == nds - 1:
        neweps = max(mstar)
        mstar = mstar * oldeps / neweps
        print oldeps, neweps
    alphavir = hst_alphavir(hstdata)
    times.append(time)
    masses.append(mstar)
    alphas.append(alphavir)

## Plotting setup
ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plow, pmid, phigh, = plt.plot(times[0], masses[0], 'k', times[1], masses[1], 'r', times[2], masses[2], 'b')
plt.legend((plow, pmid, phigh), (r"$\displaystyle {\rm N = 128}$",r"$\displaystyle {\rm N = 256}$",r"$\displaystyle {\rm N = 512}$"),prop={'size':8})
plt.axis([0,3,0,1])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*1,r"$\displaystyle(a)$")

axy = plt.subplot(1,2,2)
axy.yaxis.set_label_position("right")
plow, pmid, phigh, = plt.plot(times[0], alphas[0], 'k', times[1], alphas[1], 'r', times[2], alphas[2], 'b')
plt.axis([0,3,1.0e-1,1.0e2])
plt.yscale('log')
plt.xticks([0,1,2,3])
plt.yticks([1.0e-1,1.0,1.0e1,1.0e2])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \alpha_{\rm vir}$")
plt.text(0.05*2,10**(-1+0.9*3),r"$\displaystyle(b)$")

pp = PdfPages(plotdir + fname)
pp.savefig()
pp.close()


