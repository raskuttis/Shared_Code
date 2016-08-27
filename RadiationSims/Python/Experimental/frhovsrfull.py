from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_out import *
from hyp_field import *
from hyp_star import *
from hyp_hst import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/tigress-hsm/raskutti/tigress-hsm/Hyperion/RadParGrav/'
outdatadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
outfile = 'RadParGrav.out'
starfile = 'star'
hstfile = 'id0/RadParGrav.hst'
outhostname = 'raskutti@bellona.astro.princeton.edu'
vtkbase = 'RadParGrav'
hostname = 'raskutti@tiger.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_AllOut/'
outdatafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
Nproc = 64

outlines = read_outfile(outhostname,outdatadir + outdatafolder + outfile)
hstdata = read_hstfile(outhostname,outdatadir + outdatafolder + hstfile)
mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
rhocloud = out_rhocloud(outlines)
tff = out_tff(outlines)
time = hst_time(hstdata, tff)
stardata = read_allstars(outhostname,outdatadir + outdatafolder + starfile, time, tff)

tout = 0.6
tind = np.round(tout / 0.04).astype(int)
print tout, tind

data = AthenaDataRemote(hostname,datadir + datafolder,vtkbase,tind,Nproc)
dvol = np.prod(data.dx)
density = data.get_field('density')
xrc, yrc, zrc = field_starrecentre(data.x1zones, data.x2zones, data.x3zones, stardata, tout, time)
rad, rho, rhostd, rholow, rhohigh = rrhohist(density, xrc, yrc, zrc, rcloud, rhocloud)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.plot(rad,rho,'k')
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1,1.5])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle \rho / \rho_{\rm cl,0}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

plt.subplot(2,2,2)
plt.plot(rad,rhostd,'k')
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1,1.5])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle \rho / \rho_{\rm cl,0}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

plt.subplot(2,2,3)
plt.plot(rad,rholow,'k')
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1,1.5])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle \rho / \rho_{\rm cl,0}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

plt.subplot(2,2,4)
plt.plot(rad,rhohigh,'k')
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1,1.5])

plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.ylabel(r"$\displaystyle \rho / \rho_{\rm cl,0}$")
#plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()