from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_out import *
from hyp_hst import *
from hyp_pdf import *
from hyp_math import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'denpdfall.dat'
hostname = 'raskutti@bellona.astro.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

tff = out_tff(outlines)
msol = out_msol(outlines)
kappa = out_kappa(outlines)
dx = out_dx(outlines)
print tff
pdftime = sim_pdftime(pdflines, tff)

nlines = len(pdflines)
nts = int(np.floor(nlines / 2))
plottimes = []
taulow = []
taumid = []
tauhigh = []
for i in xrange(100, nts):
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,i-1)
    pdfcum = np.cumsum(pdfm)
    ilow = np.abs(pdfcum - 0.5).argmin()
    imid = np.abs(pdfcum - 0.9).argmin()
    ihigh = np.abs(pdfcum - 0.99).argmin()
    taulow.append(10**pdflogx[ilow] * kappa * dx / 4)
    taumid.append(10**pdflogx[imid] * kappa * dx / 4)
    tauhigh.append(10**pdflogx[ihigh] * kappa * dx / 4)
    plottimes.append(pdftime[i-1])

datafolder = 'UV_M2.0e5_R15.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
tff = out_tff(outlines)
msol = out_msol(outlines)
kappa = out_kappa(outlines)
dx = out_dx(outlines)
print tff

altpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
altpdftime = sim_pdftime(altpdflines, tff)

nlines = len(altpdflines)
nts = int(np.floor(nlines / 2))
altplottimes = []
alttaulow = []
alttaumid = []
alttauhigh = []
for i in xrange(100, nts):
    pdflogx = sim_pdfx(altpdflines)
    pdfm = sim_pdfm(altpdflines,i-1)
    pdfcum = np.cumsum(pdfm)
    ilow = np.abs(pdfcum - 0.5).argmin()
    imid = np.abs(pdfcum - 0.9).argmin()
    ihigh = np.abs(pdfcum - 0.99).argmin()
    alttaulow.append(10**pdflogx[ilow] * kappa * dx / 4)
    alttaumid.append(10**pdflogx[imid] * kappa * dx / 4)
    alttauhigh.append(10**pdflogx[ihigh] * kappa * dx / 4)
    altplottimes.append(altpdftime[i-1])

xcap = 0.05*3
ycap = 10**(-1 + 0.9 * 3)
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plow, pmid, phigh = plt.plot(plottimes,taulow,'k',plottimes,taumid,'r',plottimes,tauhigh)
#plt.legend((pfid, pnf), (r"$\displaystyle {\rm Fiducial}$",r"$\displaystyle {\rm No Feedback}$"))
plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.yscale('log')
plt.axis([0,2,1.0e-1,1.0e2])
plt.xticks([0,1,2])
plt.yticks([0.1,1.0,10.0,1.0e2])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \tau_{\rm cell}$")

plt.subplot(1,2,2)
plow, pmid, phigh = plt.plot(altplottimes,alttaulow,'k',altplottimes,alttaumid,'r',altplottimes,alttauhigh)
#plt.legend((pfid, pnf), (r"$\displaystyle {\rm Fiducial}$",r"$\displaystyle {\rm No Feedback}$"))
plt.text(xcap,ycap,r"$\displaystyle(b)$")
plt.yscale('log')
plt.axis([0,2,1.0e-1,1.0e2])
plt.xticks([0,1,2])
plt.yticks([])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle \tau_{\rm cell}$")

pp = PdfPages(plotdir + 'f7.pdf')
pp.savefig()
pp.close()


