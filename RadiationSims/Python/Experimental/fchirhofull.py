from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_out import *
from hyp_field import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/tigress-hsm/raskutti/tigress-hsm/Hyperion/RadParGrav/'
outdatadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
outfile = 'RadParGrav.out'
outhostname = 'raskutti@bellona.astro.princeton.edu'
vtkbase = 'RadParGrav'
hostname = 'raskutti@tiger.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_AllOut/'
outdatafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
Nproc = 64

outlines = read_outfile(outhostname,outdatadir + outdatafolder + outfile)
mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)

data = AthenaDataRemoteJoined(hostname,datadir + datafolder,vtkbase,15)
dvol = np.prod(data.dx)
density = data.get_field('density')
fluxrad = data.get_field('rad_flux')
xrc, yrc, zrc = field_recentre(density, data.x1zones, data.x2zones, data.x3zones)
rad, rho = rvhist(density, xrc, yrc, zrc, rcloud, 1.0)
rad, fr = rfrhist(fluxrad, xrc, yrc, zrc, rcloud, 1.0)
rad, rhofr = rfrhist(fluxrad * density, xrc, yrc, zrc, rcloud, 1.0)
chirho = rhofr / (rho * fr)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,1,1)
plt.plot(rad,chirho,'k')
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1,1.5])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()