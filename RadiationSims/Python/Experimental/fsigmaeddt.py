from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_pdf import *
import matplotlib.pyplot as plt
from hyp_math import *

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
clight = out_clight(outlines)
rcloud = out_rcloud(outlines)
msol = out_msol(outlines)
dx = out_dx(outlines)
vturb = out_vturb(outlines)
gravc = out_G(outlines)
kappa = out_kappa(outlines)

time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
eff = hst_eff(hstdata, mcloud)
psi = out_Psi(outlines)
psicgs = 2000.0

sigmaedd = surf_edd(psicgs)
sigmaeddt = sigmaedd * mstar / (1.0 + mstar)

sigmain = 1.0 * msol
eff = 0.3
vnorm = psi / (4.0 * np.pi * clight * sigmain) * mcloud / rcloud**2 * (1.0 - np.exp(-1.0 * sigmain * kappa))
vnormtwo = (psi / (4.0 * np.pi * clight * sigmain) - gravc) * mcloud / rcloud**2
eps0 = 2.0 * np.pi * clight * gravc * sigmain / (psi - 2.0 * np.pi * clight * gravc * sigmain)
t0 = eps0 / eff
plottime = np.linspace(0.0,1.0,100)
vout = vnorm * eff / 2.0 * (plottime - t0)**2 + vnorm * eps0 * (plottime - t0)
print vnormtwo, vnorm
exit()

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.plot(time,sigmaeddt)
#plt.xscale('log')
#plt.yscale('log')
plt.axis([0,3,0.0,300.0])
#plt.xticks([0.0,5.0,10.0,15.0])

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()

